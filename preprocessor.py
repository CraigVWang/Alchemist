"""The module that handles preprocessing of molecular data before alchemical calculations."""

import os
import tqdm
import re
import shutil
from typing import List, Dict, Any, Optional
from omegaconf import DictConfig
from pathlib import Path
from metadata_admin import MetadataAdmin
from converter import Converter

class Preprocessor:
    """
    结构文件预处理类
    """
    
    # 类常量 - 支持的所有格式
    SUPPORTED_EXTENSIONS = {'.pdb', '.cif', '.mol', '.mol2', '.sdf', '.xyz'}
    FILE_TYPE_MAPPING = {
        '.pdb': 'pdb',
        '.cif': 'cif',  # mmCIF 文件使用 .cif 扩展名
        '.mol': 'mol',
        '.mol2': 'mol2',
        '.sdf': 'sdf',
        '.xyz': 'xyz'
    }

    def __init__(self, *, conf: DictConfig):
        """
        初始化预处理类
        参数:
            config: Hydra配置对象
        """
        print("="*40)
        print("\n🔧 初始化预处理器...")
        # 导入配置
        self.conf = conf
        self.preprocessor_conf = self.conf.preprocessor

        # 导入目录
        self.data_dir = self.conf.data.data_dir
        self.inputs_dir = self.conf.data.input.inputs_dir
        self.raw_dir = self.conf.data.input.raw_dir
        self.outputs_dir = self.conf.data.output.outputs_dir
        self.preprocessed_dir = self.conf.data.output.preprocessed_dir
        self.ligand_dir = self.conf.data.output.ligand_dir
        self.protein_dir = self.conf.data.output.protein_dir
        self.alchemical_dir = self.conf.data.output.alchemical_dir
        self.analysis_dir = self.conf.data.output.analysis_dir
        self.visualization_dir = self.conf.data.output.visualization_dir
        self.experiments_dir = self.conf.data.experiments_dir
        self.metadata_dir = self.conf.metadata_admin.metadata_dir

        # 初始化metadata实例
        self._metadata_admin = MetadataAdmin(conf=conf)
        self._converter = Converter(conf=conf)
        self.setup_directories()
        print("📂 目录设置完成:")
        print("✅ 预处理器初始化成功!\n")
        print("="*40)

    def setup_directories(self):
        """根据配置创建必要的目录结构"""
        # 创建目录列表并遍历，依次创建目录
        dir_list = [
        self.data_dir,
        self.inputs_dir,
        self.raw_dir,
        self.outputs_dir,
        self.preprocessed_dir,
        self.ligand_dir,
        self.protein_dir,
        self.alchemical_dir,
        self.analysis_dir,
        self.visualization_dir,
        self.experiments_dir,
        self.metadata_dir,
        ]

        for d in dir_list:
            dir = Path(d)
            if not dir.exists():
                dir.mkdir(parents=True,exist_ok=True)
                print(f"📁 创建原始数据目录: {dir}")

        print("📂 目录设置完成:")
        print(f"   - 数据目录: {self.data_dir}")
        print(f"   - 输入处理目录: {self.inputs_dir}")
        print(f"   - 原始数据目录: {self.raw_dir}")
        print(f"   - 输出数据文件: {self.outputs_dir}")
        print(f"   - 预处理数据文件: {self.preprocessed_dir}")
        print(f"   - 配体分子文件: {self.ligand_dir}")
        print(f"   - 蛋白质文件: {self.protein_dir}")
        print(f"   - 炼金术数据文件: {self.alchemical_dir}")
        print(f"   - 分析数据文件: {self.analysis_dir}")
        print(f"   - 可视化数据文件: {self.visualization_dir}")
        print(f"   - 实验数据文件: {self.experiments_dir}")
        print(f"   - 元数据文件: {self.metadata_dir}")

    def extract_pdb_id(self, filename: str) -> str:
        """
        从文件名中提取可能的PDB ID
        
        PDB ID通常是4个字符的代码，第一个是数字1-9，后面三个是字母或数字
        
        参数:
            filename: 文件名
            
        返回:
            PDB ID字符串或'NAN'
        """
        pdb_pattern = r'[1-9][a-z0-9]{3}'
        matches = re.findall(pdb_pattern, filename.lower())
        
        for match in matches:
            if len(match) == 4 and re.match(r'^[1-9a-z][a-z0-9]{3}$', match):
                return match.upper()
        
        return 'NAN'
    
    def get_file_type(self, filename: str) -> str:
        """
        根据文件扩展名确定文件类型
        
        参数:
            filename: 文件名
            
        返回:
            文件类型字符串
        """
        ext = Path(filename).suffix.lower()
        return self.FILE_TYPE_MAPPING.get(ext, 'unknown')

    
    def scan_directory(self, existing_metadata: List[Dict[str, str]]) -> List[Dict[str, str]]:
        """
        扫描目录中的结构文件，与现有元数据合并
        
        参数:
            existing_metadata: 现有元数据列表
            
        返回:
            包含文件信息的字典列表
        """
        root_dir = self.raw_dir
        
        # 创建现有分子的查找字典
        existing_molecules = {item['name']: item for item in existing_metadata}
        new_data = []

        print("="*40)
        print(f"🔍 处理所有支持格式的文件")
        supported_extensions = self.SUPPORTED_EXTENSIONS
              
        print(f"📁 开始扫描目录: {root_dir}")
        
        # 收集所有文件
        all_files = []
        for root, dirs, files in os.walk(root_dir):
            for file in files:
                all_files.append((root, file))
        
        # 使用进度条显示扫描进度
        file_iterator = tqdm(all_files, desc="📂 扫描文件")
        
        # 过滤支持的文件格式
        for root, file in file_iterator:
            file_path = os.path.join(root, file)
            file_ext = Path(file).suffix.lower()
            
            if file_ext in supported_extensions:
                mol_name = Path(file).stem
                
                # 如果分子已经在元数据中，跳过扫描（只更新新文件）
                if mol_name in existing_molecules:
                    continue
                
                pdb_id = self.extract_pdb_id(mol_name)
                file_type = self.get_file_type(file)
                
                # 计算相对于原始目录的相对路径
                relative_path = os.path.relpath(root, root_dir)
                
                # 创建新分子的默认数据
                new_data.append({
                    'name': mol_name,
                    'filename': file,
                    'original_file_path': file_path,
                    'relative_path': relative_path,
                    'pdb_id': pdb_id,
                    'original_file_type': file_type,
                    
                    # 初始化所有列
                    'preprocessed_file_path': '',
                    'preprocessed_file_type': '',
                    'prepared_system_path': '',
                    'alchemical_result_path': '',
                    'analysis_result_path': '',
                    
                    'preprocess_success': 'False',
                    'preparation_success': 'False',
                    'alchemical_success': 'False',
                    'analysis_success': 'False',
                    'finish_success': 'False',
                    
                    'preprocess_timestamp': '',
                    'preparation_timestamp': '',
                    'alchemical_timestamp': '',
                    'analysis_timestamp': '',
                    
                    'ligand_atom_count': '0',
                    'free_energy_value': '0.0',
                    'free_energy_error': '0.0',
                    'processing_notes': ''
                })
        
        # 合并现有数据和新数据
        merged_data = existing_metadata + new_data
        
        print(f"✅ 找到 {len(merged_data)} 个分子（{len(existing_molecules)} 个现有 + {len(new_data)} 个新）")
        print("="*40)
        return merged_data
    
    def process_molecule_file(self, input_path: str, mol_name: str, file_type: str, relative_path: str) -> Optional[str]:
        """
        处理分子文件 - 复制或转换格式，保持目录结构
        
        参数:
            input_path: 输入文件路径
            mol_name: 分子名称
            file_type: 文件类型
            relative_path: 相对于原始目录的路径
            
        返回:
            处理后的文件路径，如果失败返回None
        """
        # 构建输出目录路径，保持原始目录结构
                
        mol_output_dir = Path(self.ligand_dir) / relative_path
        mol_output_dir.mkdir(parents=True, exist_ok=True)
        pdb_output_dir = Path(self.protein_dir) / relative_path
        pdb_output_dir.mkdir(parents=True, exist_ok=True)
        
        # 确定输出文件路径和格式
        if file_type in ['pdb', 'cif',]:
            # 对于这些格式，保持原格式，直接复制
            pdb_output_path = pdb_output_dir / Path(input_path).name
            try:
                shutil.copy2(input_path, pdb_output_path)
                return str(pdb_output_path), None
            except Exception as e:
                print(f"❌ 复制失败 {input_path}: {e}")
                return None
        
        if file_type in ['mol', 'mol2', 'sdf']:
            # 对于这些格式，保持原格式，直接复制
            mol_output_path = mol_output_dir / Path(input_path).name
            try:
                shutil.copy2(input_path, mol_output_path)
                return str(mol_output_path), None
            except Exception as e:
                print(f"❌ 复制失败 {input_path}: {e}")
                return None
            
        elif file_type == 'xyz':
            # XYZ 文件转换为 MOL2 格式
            mol_path, pdb_path = self._converter.convert(input_path, mol_name, file_type, relative_path)
            return mol_path, pdb_path
        
        else:
            print(f"⚠️ 跳过不支持的文件格式: {input_path}")
            return None
        
    
    def batch_process_files(self, metadata: List[Dict[str, str]], test_single: bool = False) -> List[Dict[str, str]]:
        """
        批量处理文件，保持目录结构
        
        参数:
            metadata: 元数据列表
            test_single: 是否只测试单个样本
            
        返回:
            更新后的元数据列表
        """
        successful_processing = 0
        total_to_process = len(metadata)
        
        # 如果测试单个样本，只处理第一个
        if test_single and metadata:
            metadata = [metadata[0]]
            total_to_process = 1
            print(f"🧪 测试模式：只处理第一个分子: {metadata[0]['name']}")
        
        # 使用进度条显示处理进度
        data_iterator = tqdm(metadata, desc="🔄 处理文件")
        
        for item in data_iterator:
            mol_name = item['name']
            input_file = item['original_file_path']
            file_type = item['original_file_type']
            relative_path = item['relative_path']
            
            # 如果已经处理成功，跳过
            if item.get('preprocess_success', 'False').lower() == 'true':
                data_iterator.set_postfix_str(f"跳过: {successful_processing}/{total_to_process}")
                continue
            
            # 处理文件，传递相对路径
            output_file = self.process_molecule_file(input_file, mol_name, file_type, relative_path)
                           
            if output_file:
                for i in range(2):
                    # 确定输出文件类型
                    if output_file[i] is not None:
                        output_path = Path(output_file[i])
                        output_file_type = output_path.suffix[1:]  # 去掉点号
                        
                        # 更新元数据
                        self._metadata_admin.update_molecule_status(
                            metadata=metadata,
                            molecule_name=mol_name,
                            stage='preprocess',
                            success=True,
                            additional_info={
                                'preprocessed_file_path': output_file,
                                'preprocessed_file_type': output_file_type
                            }
                        )
                        
                        successful_processing += 1
                        # 更新进度条描述
                        data_iterator.set_postfix_str(f"成功: {successful_processing}/{total_to_process}")
                    else:
                        continue

            else:
                print(f"❌ 处理失败: {input_file}")
                # 更新状态为失败
                self._metadata_admin.update_molecule_status(
                    metadata=metadata,
                    molecule_name=mol_name,
                    stage='preprocess',
                    success=False
                )
        
        print(f"📊 成功处理 {successful_processing}/{total_to_process} 个文件")
        return metadata

    def run(self, test_single: bool = False) -> Dict[str, Any]:
        """
        主要的预处理流程
        
        参数:
            test_single: 是否只测试单个样本
            
        返回:
            处理结果的字典，包含统计信息和文件路径
        """
        print("=" * 40)
        print("🚀 开始结构文件预处理流程")
        if test_single:
            print("🧪 测试模式：只处理单个样本")
        print("=" * 40)
        
        # 初始化元数据文件
        metadata = self._metadata_admin.initialize_metadata_file()
        
        # 扫描目录并更新元数据
        metadata = self.scan_directory(metadata)
        
        if not metadata:
            print("❌ 未找到任何结构文件")
            return {'success': False, 'message': '未找到任何结构文件'}
        
        # 执行文件处理
        print("\n🔄 开始批量处理文件...")
        metadata = self.batch_process_files(metadata, test_single)
        
        # 保存元数据
        self._metadata_admin.save_metadata(metadata)
        
        # 生成统计信息
        stats = self._metadata_admin.generate_statistics(metadata)
        self._metadata_admin.print_statistics(stats)
        
        result = {
            'success': True,
            'total_molecules': len(metadata),
            'preprocess_success': stats['preprocess_success'],
            'preprocess_success_rate': stats['preprocess_success']/len(metadata) if len(metadata) > 0 else 0,
            'metadata_file': str(self.metadata_file),
            'metadata': metadata  # 返回元数据供后续阶段使用
        }
        
        print("=" * 40)
        print("🎉 结构文件预处理完成!")
        print("=" * 40)
        return result    
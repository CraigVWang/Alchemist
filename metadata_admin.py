"""The module that handles metadata administration."""

import csv
from pathlib import Path
from omegaconf import DictConfig
from datetime import datetime

class MetadataAdmin:
    """
    本模块负责元数据的管理和维护。包括元数据的创建、更新、删除和查询等功能。
    """
    def __init__(self, *, conf: DictConfig):
        """
        初始化MetadataAdmin类。
        """
        self.conf = conf
        self.metadata_admin_conf = self.conf.metadata_admin
        self.metadata_file = Path(self.metadata_admin_conf.metadata_dir) / "metadata.csv"
        self.metadata = []

    def load_metadata(self) -> list:
        """
        加载元数据文件
        
        返回:
            元数据列表
        """
        if not self.metadata_file.exists():
            print(f"❌ 元数据文件不存在: {self.metadata_file}")
            return []
        
        try:
            with open(self.metadata_file, 'r', newline='', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                self.metadata = list(reader)
            print(f"📖 加载元数据: {len(self.metadata)} 个分子")
            return self.metadata
        except Exception as e:
            print(f"❌ 加载元数据失败: {e}")
            return []
    
    def save_metadata(self):
        """
        保存元数据到文件
        """
        if not self.metadata:
            print(f"💾 元数据是空文件")
            
        try:
            # 获取所有可能的列
            all_columns = set()
            for item in self.metadata:
                all_columns.update(item.keys())
            
            # 写入文件
            with open(self.metadata_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=sorted(all_columns))
                writer.writeheader()
                for row in self.metadata:
                    writer.writerow(row)
            
            print(f"💾 元数据已保存到: {self.metadata_file}")
        except Exception as e:
            print(f"❌ 保存元数据失败: {e}")

    def initialize_metadata_file(self):
        """
        初始化元数据文件，如果不存在则创建，如果存在则读取
        """
        if not self.metadata_file.exists():
            # 创建新文件
            columns = [
                'name',                     # 分子名称
                'filename',                 # 原始文件名
                'original_file_path',       # 原始文件路径
                'relative_path',            # 相对路径
                'pdb_id',                   # PDB ID
                'original_file_type',       # 原始文件类型
                'preprocessed_file_path',   # 预处理后文件路径
                'preprocessed_file_type',   # 预处理后文件类型
                'prepared_system_path',     # 准备系统路径
                'alchemical_result_path',   # 炼金术结果路径
                'analysis_result_path',     # 分析结果路径
                
                # 状态列
                'preprocess_success',   # 预处理是否成功
                'preparation_success',   # 最小化是否成功
                'alchemical_success',  # 炼金术是否成功
                'analysis_success',    # 分析是否成功
                'finish_success',      # 全部完成是否成功
                
                # 时间戳
                'preprocess_timestamp',     # 预处理时间
                'preparation_timestamp',    # 系统准备时间
                'alchemical_timestamp',     # 炼金术时间
                'analysis_timestamp',       # 分析时间
                
                # 统计信息
                'ligand_atom_count',        # 配体原子数
                'free_energy_value',        # 自由能值
                'free_energy_error',        # 自由能误差
                'processing_notes'          # 处理备注
            ]
            
            with open(self.metadata_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=columns)
                writer.writeheader()
            
            print(f"📄 创建新元数据文件: {self.metadata_file}")
            return []
        else:
            # 读取现有文件
            print(f"📖 读取现有元数据文件: {self.metadata_file}")
            with open(self.metadata_file, 'r', newline='', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                return list(reader)
            
    def generate_statistics(self, metadata: List[Dict[str, str]]) -> Dict[str, Any]:
        """
        生成文件统计信息
        
        参数:
            metadata: 元数据列表
            
        返回:
            统计信息字典
        """
        stats = {
            'total_molecules': len(metadata),
            'file_types': {},
            'pdb_ids_count': 0,
            'processed_success': 0,
            'preparation_success': 0,
            'alchemical_success': 0,
            'analysis_success': 0,
            'finish_success': 0
        }
        
        for item in metadata:
            # 文件类型统计
            file_type = item['original_file_type']
            stats['file_types'][file_type] = stats['file_types'].get(file_type, 0) + 1
            
            # PDB ID统计
            if item.get('pdb_id', 'NAN') != 'NAN':
                stats['pdb_ids_count'] += 1
            
            # 状态统计
            if item.get('processed_successfully', 'False').lower() == 'true':
                stats['processed_success'] += 1
            
            if item.get('minimized_successfully', 'False').lower() == 'true':
                stats['preparation_success'] += 1
            
            if item.get('alchemical_successfully', 'False').lower() == 'true':
                stats['alchemical_success'] += 1
            
            if item.get('analysis_successfully', 'False').lower() == 'true':
                stats['analysis_success'] += 1
            
            if item.get('finish_successfully', 'False').lower() == 'true':
                stats['finish_success'] += 1
        
        return stats
    def update_molecule_status(self, 
                              molecule_name: str,
                              stage: str,
                              success: bool = True,
                              additional_info: dict = None):
        """
        更新分子的状态信息
        
        参数:
            molecule_name: 分子名称
            stage: 阶段名称 ('preprocess', 'preparation', 'alchemical', 'analysis')
            success: 该阶段是否成功
            additional_info: 额外的信息字典
        """
        if not self.metadata:
            return
            
        # 查找分子
        for item in self.metadata:
            if item['name'] == molecule_name:
                # 更新时间戳
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                
                # 更新对应阶段的状态
                if stage == 'preprocess':
                    item['preprocess_success'] = str(success)
                    item['preprocess_timestamp'] = timestamp
                
                elif stage == 'preparation':
                    item['preparation_success'] = str(success)
                    item['preparation_timestamp'] = timestamp
                
                elif stage == 'alchemical':
                    item['alchemical_success'] = str(success)
                    item['alchemical_timestamp'] = timestamp
                
                elif stage == 'analysis':
                    item['analysis_success'] = str(success)
                    item['analysis_timestamp'] = timestamp
                    
                    # 计算finish_successfully
                    conditions = [
                        item.get('preprocess_success', 'False').lower() == 'true',
                        item.get('preparation_success', 'False').lower() == 'true',
                        item.get('alchemical_success', 'False').lower() == 'true',
                        success  # 当前的分析阶段是否成功
                    ]
                    finish_success = all(conditions)
                    item['finish_successfully'] = str(finish_success)
                
                # 更新额外信息
                if additional_info:
                    for key, value in additional_info.items():
                        item[key] = str(value) if value is not None else ''
                
                break
    
    def get_statistics(self) -> dict:
        """
        获取处理统计信息
        
        返回:
            统计信息字典
        """
        if not self.metadata:
            return {
                'total_molecules': 0,
                'preprocess_success': 0,
                'preparation_success': 0,
                'alchemical_success': 0,
                'analysis_success': 0,
                'finish_success': 0
            }
        
        stats = {
            'total_molecules': len(self.metadata),
            'preprocess_success': 0,
            'preparation_success': 0,
            'alchemical_success': 0,
            'analysis_success': 0,
            'finish_success': 0
        }
        
        for item in self.metadata:
            if item.get('preprocess_success', 'False').lower() == 'true':
                stats['preprocess_success'] += 1
            
            if item.get('preparation_success', 'False').lower() == 'true':
                stats['preparation_success'] += 1
            
            if item.get('alchemical_success', 'False').lower() == 'true':
                stats['alchemical_success'] += 1
            
            if item.get('analysis_success', 'False').lower() == 'true':
                stats['analysis_success'] += 1
            
            if item.get('finish_success', 'False').lower() == 'true':
                stats['finish_success'] += 1
        
        # 计算成功率
        for key in ['preprocess', 'preparation', 'alchemical', 'analysis', 'finish']:
            total_key = 'total_molecules'
            success_key = f'{key}_success'
            if success_key in stats and stats[total_key] > 0:
                stats[f'{key}_success_rate'] = stats[success_key] / stats[total_key]
        
        return stats
    
    def print_statistics(self):
        """打印处理统计信息"""
        stats = self.get_statistics()
        
        print("\n📊 处理统计信息:")
        print("=" * 40)
        print(f"总分子数: {stats['total_molecules']}")
        print(f"预处理成功: {stats['preprocess_success']} ({stats.get('preprocess_success_rate', 0)*100:.1f}%)")
        print(f"系统准备成功: {stats['preparation_success']} ({stats.get('preparation_success_rate', 0)*100:.1f}%)")
        print(f"炼金术成功: {stats['alchemical_success']} ({stats.get('alchemical_success_rate', 0)*100:.1f}%)")
        print(f"分析成功: {stats['analysis_success']} ({stats.get('analysis_success_rate', 0)*100:.1f}%)")
        print(f"完成全部流程: {stats['finish_success']} ({stats.get('finish_success_rate', 0)*100:.1f}%)")
        print("=" * 40)
              


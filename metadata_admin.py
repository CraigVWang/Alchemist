"""The module that handles metadata administration."""

import csv
from pathlib import Path
from omegaconf import DictConfig

class MetadataAdmin:
    """
    本模块负责元数据的管理和维护。包括元数据的创建、更新、删除和查询等功能。
    """
    def __init__(self, *, conf: DictConfig):
        """
        初始化MetadataAdmin类。
        """
        self.metadata_admin_conf = conf
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

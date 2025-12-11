"""
The module that create the pipeline of the whole workflow.
Version:1.0
Author:CraigVWang
"""

import hydra
from omegaconf import DictConfig, OmegaConf

# Import internal modules
from metadata_admin import MetadataAdmin
from converter import Converter
from preprocessor import Preprocessor
from provider import Provider
from alchemist import Alchemist
from analyzer import Analyzer


class Commander:
    def __init__(self, *, conf: DictConfig):
        """
        初始化Commander类，设置各个模块的配置。
        参数:
            conf (DictConfig): 包含各个模块配置的字典配置对象。
        """
        self.conf = conf
        self.metadata_admin_conf = self.conf.metadata_admin
        self.converter_conf = self.conf.converter
        self.preprocessor_conf = self.conf.preprocessor
        self.provider_conf = self.conf.provider
        self.alchemist_conf = self.conf.alchemist
        self.analyzer_conf = self.conf.analyzer

        # 各个模块的实例设置为None，后续调用时初始化
        self.preprocessor = None
        self.provider = None
        self.alchemist = None
        self.analyzer = None
        
        # 初始化并统一元数据
        self.metadata = []
        
    def run_preprocessor(self, test_single: bool = False):
        """运行预处理模块"""
        print("="*40)
        print("🚀 开始数据预处理阶段")
        if test_single:
            print("🧪 测试模式：只处理单个样本")
        print("="*40)

        if self.preprocessor is None:
            self.preprocessor = Preprocessor(conf=self.conf)
        
        result = self.preprocessor.run(test_single)
        self.metadata = result['metadata']

    def run_provider(self, test_single: bool = False):
        """运行系统准备模块"""
        print("="*40)
        print("🔄 运行系统准备模块...")
        print("="*40)

        if self.provider is None:
            self.provider = Provider(conf=self.conf)
        
        self.provider.run_preparation(self.metadata, test_single)
        
    def run_alchemist(self, test_single: bool = False):
        """运行炼金术模块"""
        print("="*40)
        print("🔄 运行炼金术模块...")
        print("="*40)

        if self.alchemist is None:
            self.alchemist = Alchemist(conf=self.conf)
        

    def run_analyzer(self, test_single: bool = False):
        """运行分析模块"""
        print("="*40)
        print("🔄 运行分析模块...")
        print("="*40)

        if self.analyzer is None:
            self.analyzer = Analyzer(conf=self.conf)
        

    def run_preprocessing_only(self):
        """仅运行预处理步骤"""
        self.run_preprocessor()

    def run_full_pipeline(self):
        """运行完整的工作流程"""
        print("="*40)
        print("🚀 开始数据预处理、系统准备和炼金术阶段")
        print("="*40)
        
        self.run_preprocessor()
        self.run_provider()
        self.run_alchemist()
        self.run_analyzer()

    def run_single_test(self,test_single: bool = True):
        """运行单一测试"""
        print("="*40)
        print("🔄 运行单一测试...")
        print("="*40)
        
        self.run_preprocessor(test_single)
        self.run_provider(test_single)
        self.run_alchemist(test_single)
        self.run_analyzer(test_single)
        

@hydra.main(version_base=None, config_path="./config", config_name="base")
def main(cfg: DictConfig):
    """主函数"""
    print("⚙️ 实验配置:")
    print(OmegaConf.to_yaml(cfg))
    
    # 创建实验实例
    experiment = Commander(cfg)
    
    # 根据配置选择运行模式
    mode = cfg.commander.get('mode', 'full')
    
    if mode == 'preprocess_only':
        experiment.run_preprocessing_only()
    elif mode == 'full':
        experiment.run_full_pipeline()
    elif mode == 'test_single':
        experiment.run_single_test()
    else:
        print(f"❌ 未知的运行模式: {mode}")
        print("可用模式: preprocess_only, full, test_single")

if __name__ == "__main__":
    main()
"""The module that create the pipeline of the whole workflow."""

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

        # 初始化各个模块的实例
        self.metadata_admin = MetadataAdmin(conf=self.metadata_admin_conf)
        self.converter = Converter(conf=self.converter_conf)
        self.preprocessor = Preprocessor(conf=self.preprocessor_conf)
        self.provider = Provider(conf=self.provider_conf)
        self.alchemist = Alchemist(conf=self.alchemist_conf)
        self.analyzer = Analyzer(conf=self.analyzer_conf)
        
    def run_converter(self):
        """运行转换模块"""
        print("🔄 运行转换模块...")
        # 在这里添加转换逻辑
        pass

    def run_preprocessor(self):
        """运行预处理模块"""
        print("🔄 运行预处理模块...")
        # 在这里添加预处理逻辑
        pass

    def run_provider(self):
        """运行系统准备模块"""
        print("🔄 运行系统准备模块...")
        # 在这里添加提供者逻辑
        pass

    def run_alchemist(self):
        """运行炼金术模块"""
        print("🔄 运行炼金术模块...")
        # 在这里添加Alchemist逻辑
        pass

    def run_analyzer(self):
        """运行分析模块"""
        print("🔄 运行分析模块...")
        # 在这里添加分析逻辑
        pass

    def run_preprocessing_only(self):
        """仅运行预处理步骤"""
        self.run_converter()
        self.run_preprocessor()

    def run_full_pipeline(self):
        """运行完整的工作流程"""
        self.run_converter()
        self.run_preprocessor()
        self.run_provider()
        self.run_alchemist()
        self.run_analyzer()

    def run_single_test(self):
        """运行单一测试"""
        print("🔄 运行单一测试...")
        # 在这里添加单一测试逻辑
        pass


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
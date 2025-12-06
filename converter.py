"""The module that handles conversion of the molecular and protein data formats."""

from omegaconf import DictConfig

class Converter:
    """
    本模块负责分子和蛋白质数据格式的转换。将所有的输入数据转换为统一的格式，以便后续处理。预计支持的格式包括PDB、MMCIF、MOL2、SDF、XYZ。
    """
    def __init__(self, *, conf: DictConfig):
        """
        初始化Converter类，设置转换模块的配置。
        参数:
            conf (DictConfig): 包含转换模块配置的字典配置对象。
        """
        self.converter_conf = conf

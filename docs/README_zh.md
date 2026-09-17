# PACE 中文文档

请从仓库的[中文使用说明](../README.zh-CN.md)进入。

- [命令行与三模式真实数据参数](cli.md)
- [当前总公式](FORMULA.md)与[模型说明](model.md)
- [连续 eta 校准](eta_calibration.md)
- [输入准备](input_preparation.md)、[参数](parameters.md)和[数据字典](data_dictionary.md)
- [训练](training.md)、[比较与基准](comparison.md)、[验证记录](validation.md)

全部入口对应同一个已安装的 PACE 实现。默认无适用功能校准时 eta=0；
有适用训练/校准标签时估计 [0,1] 中的连续值。测试集不参与拟合。

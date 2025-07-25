### 一、代码阅读
原代码仅有一种单元————三维两节点线性单元

识别节点类型的原理是通过`ElementGroup`类内结合`ElementTypes`实现的，其中关键在于
```C++
void CElementGroup::CalculateMemberSize()

void CElementGroup::AllocateElements(std::size_t size)

void CElementGroup::AllocateMaterials(std::size_t size)
```
这三个函数中的`switch`语句

而具体读取单元数据的是`CBar`类中继承的`bool CBar::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)`函数

关键在于完成`CT3`



### 二、修改目标
#### 增添三维三节点单元`CT3`
正在进行，已添加头文件`T3.h`和`T3.cpp`。

已完成`T3.h`。

#### 在`CMaterial`中添加三维三节点单元材料类`CT3Material`
注意到`Matrial.h`中也需要为具体材料添加代码。
已完成，相较于基础类添加了`Thickness`厚度作为材料属性，此外，还添加了泊松比`nu`和平面应力指示变量`plane_stress`。

#### 在`ElementGroup`类中完善`switch`语句

#### 在`Outputter`的`switch`语句完善关于应力输出

#### 分片实验表
|节点号$I$|$x_I$|$y_I$|$u_I$|$v_I$|$F_{xI}$|$F_{yI}$|
|--|--|--|--|--|--|--|
|1|0|0|0|0|-10|0|
|2|2.5|0|0.025|0|15|0|
|3|2.5|3|0.025|-0.009|10|0|
|4|0|2|0|-0.006|-15|0|
|5|1|1.2|0.01|-0.0036|0|0|

### 三、TODO
- [x] 完成`CT3Material`的`ComputeElasticityMatrix()`方法
- [x] 完成`CT3Material`的`ElementStress()`方法
- [x] 在`OutPutter`类的`OutputElementStress()`方法中的`switch`语句中添加`T3`单元
- [x] 完成`T3`的`GenerateLocationMatrix()`方法（好像源代码已经完成了
- [x] 完成`ElementGroup`类的三个`switch`
- [x] 针对`T3`实现`GenerateLocationMatrix()`。要求只有6个单元
- [x] 明确2D和3D的问题
- [x] 尝试编译
- [x] 编辑测试文件并测试（例题4-4的三角形单元版本）
- [x] 利用Abaquas验证，位移计算结果一致，但是输出应力时`Segmentation fault (core dumped)`
- [x] 完善`Outputter`类中的的`OutputT3Elements`
- [x] debug，使其能正确输出应力
- [x] 后处理绘图输出位移云图
- [x] 收敛率分析（通过4个不同的h计算误差范数拟合双对数曲线斜率），参考例题4.12
- [x] 原有代码不支持添加非零的位移边界条件，导致难以进行分片实验的A和B；需要添加非零边界条件支持
- [x] 增加节点力到输出
- [x] 处理无节点力输入、无自由点等特殊情况
- [x] 分片实验（基于例题4-5）
- [x] 完成课程设计报告
- [ ] 整理代码库文件结构
### 代码阅读
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



### 修改目标
#### 增添三维三节点单元
正在进行，已添加头文件`T3.h`和`T3.cpp`。

#### 在原有代码中添加三维三节点单元
注意到`Matrial.h`中也需要为具体材料添加代码。
已完成，相较于基础类添加了`Thickness`厚度作为材料属性。
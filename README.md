# STAPpp with T3 element

## Project Overview

This is an enhanced version of the STAPpp finite element analysis program. Building upon the [original STAPpp](https://github.com/Maxwowq/STAPpp) framework, this project introduces **T3 triangular element support** and implements several significant functionality extensions.

## 🚀 Key Improvements & New Features

### 1. T3 Triangular Element Implementation
- **New CT3 Element Class**: Complete implementation of 3-node triangular elements
- **CT3Material Class**: Support for plane stress/plane strain analysis
- **Element Stiffness Matrix Calculation**: Efficient computation based on B-matrix and elasticity matrix
- **Stress Calculation**: Accurate element stress analysis and output

### 2. Non-Zero Essential Boundary Condition Support
- **Enhanced Boundary Condition Handling**: Support for non-zero displacement boundary conditions
- **gbcode Numbering System**: Optimized degree-of-freedom numbering strategy
- **GlobalLocationMatrix**: Global location matrix generation
- **EssentialBoundary Class**: Specialized handling of essential boundary conditions

### 3. Nodal Force Computation & Output
- **Nodal Reaction Force Calculation**: Node force analysis based on stiffness matrix and displacements
- **ComputeNodalForce Method**: Efficient nodal force calculation algorithm
- **Complete Force Balance Verification**: Support for force balance checking in patch tests

### 4. Tecplot Visualization Support
- **CTecplot Class**: Dedicated Tecplot data output class
- **Displacement Contour Data**: Automatic generation of Tecplot format visualization data
- **Deformed Mesh Display**: Support for deformation visualization with magnification factors

## 📁 Project Structure

```
STAPpp/
├── src/                    # Source code directory
├── data/                   # Test data files
│   ├── patchtest/         # Patch test data
│   ├── convergence/       # Convergence analysis data
│   └── validation/        # Validation test cases
├── docs/                   # Auto-generated documentation
├── report/                 # Project reports and analysis
└── make/                   # Build configuration
```

## 🛠️ Usage
```bash
./stap++ input_file.dat
```
***Note:​​ The stap++ binary in the make directory is compiled for Linux and compatibility with other systems is not guaranteed.***

### Input File Format
```
Title
NUMNP   NUMEG   NLCASE  MODEX
Node_Number   bcode_1   bcode_2   bcode_3   X   Y   Z
Num_NoneZeroEssentialBoundary
Node_Number   Direction_Number   Displacement
LoadCase_Number   Num_Loads
Node_Number   Direction_Number   Magnitude
ElementType   Num_Element   Num_Material
Material_Number   E   nu   t   PlaneStress
Element_Number   Node1   Node2   Node3
```

## 📚 References
- Zhang, Xiong. Fundamentals of Finite Element Method. Beijing: Higher Education Press, 2023.
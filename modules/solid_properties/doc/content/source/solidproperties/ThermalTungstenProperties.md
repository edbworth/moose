# ThermalTungstenProperties

!syntax description /SolidProperties/ThermalTungstenProperties

## Description

This userobject provides thermal properties for tungsten as a function of temperature using 
correlations from [!cite](milner2024space).

!include solid_properties_units.md

## Thermal Conductivity

For temperature range: 1 ≤ T < 55 K, the thermal conductivity is defined as

\begin{equation}
k(T) =
\frac{A_0 \cdot \left( \frac{T}{1000} \right)^{0.874}}{1 + A_1 \cdot \left( \frac{T}{1000} \right) + A_2 \cdot \left( \frac{T}{1000} \right)^2 + A_3 \cdot \left( \frac{T}{1000} \right)^3}
\end{equation}

where $k(T)$ is the thermal conductivity [W/(m·K)] and $T$ is temperature [K].

!table caption=Thermal conductivity parameters for T < 55 K
| Constant | Value      |
|----------|------------|
| $A_0$    | 7.348×10⁵  |
| $A_1$    | 2.544×10¹  |
| $A_2$    | -8.304×10³ |
| $A_3$    | 1.180×10⁶  |

For temperature range 55 ≤ T ≤ 3653 K, the thermal conductivity is defined as

\begin{equation}
k(T) = \frac{B_0 + B_1 \cdot \left( \frac{T}{1000} \right) + B_2 \cdot \left( \frac{T}{1000} \right)^2 + B_3 \cdot \left( \frac{T}{1000} \right)^3}{C_0 + C_1 \cdot \left( \frac{T}{1000} \right) + \left( \frac{T}{1000} \right)^2}
\end{equation}

!table caption=Thermal conductivity parameters for T ≥ 55 K
| Constant | Value       |
|----------|-------------|
| $B_0$    | -3.679      |
| $B_1$    | 1.181×10²   |
| $B_2$    | 5.879×10¹   |
| $B_3$    | 2.867       |
| $C_0$    | -2.052×10⁻² |
| $C_1$    | 4.741×10⁻¹  |

## Specific Heat

For temperature range: 11 ≤ T ≤ 293 K, the isobaric specific heat capacity is described as:

\begin{equation}
C_p(T) =\frac{A_0 \cdot \left( \frac{T}{1000} \right)^{3.03}}{1 + A_1 \cdot \left( \frac{T}{1000} \right) + A_2 \cdot \left( \frac{T}{1000} \right)^2 + A_3 \cdot \left( \frac{T}{1000} \right)^3}
\end{equation}

where $C_p(T)$ is the specific heat [J/(kg·K)] and $T$ is temperature [K].

**Note on units**: The NASA source [!cite](milner2024space) reports specific heat in J/(g·K). 
These correlations have been converted to SI units [J/(kg·K)] by multiplying by 1000, consistent 
with the solid_properties module standard.

!table caption=Specific heat parameters for T ≤ 293 K
| Constant | Value (J/g·K) |
|----------|---------------|
| $A_0$    | 3.103×10²     |
| $A_1$    | -8.815        |
| $A_2$    | 1.295×10²     |
| $A_3$    | 1.874×10³     |

For temperature range: 293 < T ≤ 3700 K, the specific heat is defined as

\begin{equation}
C_p(T) = B_0 + B_1 \cdot \left( \frac{T}{1000} \right) + B_2 \cdot \left( \frac{T}{1000} \right)^2 + B_3 \cdot \left( \frac{T}{1000} \right)^3 + \frac{B_{-2}}{\left( \frac{T}{1000} \right)^2}
\end{equation}

!table caption=Specific heat parameters for T > 293 K
| Constant   | Value (J/g·K) |
|------------|---------------|
| $B_0$      | 1.301×10⁻¹    |
| $B_1$      | 2.225×10⁻²    |
| $B_2$      | -7.224×10⁻³   |
| $B_3$      | 3.539×10⁻³    |
| $B_{-2}$   | -3.061×10⁻⁴   |

## Density

For the temperature range: 5 ≤ T ≤ 3600 K, the density is defined as

\begin{equation}
\rho(T) = \frac{\rho_0}{\left(1 + \frac{\Delta L / L_0(T)}{100}\right)^3}
\end{equation}

where $\rho(T)$ is the density [kg/m³], $\rho_0 = 19250$ kg/m³ is the reference density, 
and the thermal expansion $\Delta L / L_0(T)$ [%] is given by:

\begin{equation}
\frac{\Delta L(T)}{L_0} = A_0 + A_1 \cdot \left( \frac{T}{1000} \right) + A_2 \cdot \left( \frac{T}{1000} \right)^2 + A_3 \cdot \left( \frac{T}{1000} \right)^3
\end{equation}

!table caption=Thermal expansion parameters for T ≤ 294 K
| Constant | Value      |
|----------|------------|
| $A_0$    | -8.529×10⁻² |
| $A_1$    | -9.915×10⁻² |
| $A_2$    | 2.257      |
| $A_3$    | -3.157     |

!table caption=Thermal expansion parameters for T > 294 K
| Constant | Value      |
|----------|------------|
| $A_0$    | -1.400×10⁻¹ |
| $A_1$    | 4.869×10⁻¹  |
| $A_2$    | -3.056×10⁻² |
| $A_3$    | 2.234×10⁻²  |

## Range of Validity

- Thermal conductivity: 1 K ≤ T ≤ 3653 K (-272°C to 3380°C)
- Specific heat: 11 K ≤ T ≤ 3700 K (-262°C to 3427°C)
- Density: 5 K ≤ T ≤ 3600 K (-268°C to 3327°C)

## Usage

This userobject is typically used with [ThermalSolidPropertiesMaterial](/ThermalSolidPropertiesMaterial.md) 
to evaluate properties at quadrature points:

```
[SolidProperties]
  [tungsten]
    type = ThermalTungstenProperties
  []
[]

[Materials]
  [tungsten_mat]
    type = ThermalSolidPropertiesMaterial
    sp = tungsten
    temperature = T
  []
[]
```

!syntax parameters /SolidProperties/ThermalTungstenProperties

!syntax inputs /SolidProperties/ThermalTungstenProperties

!syntax children /SolidProperties/ThermalTungstenProperties

!bibtex bibliography

The project is still in development.

# RFSolver

A program that calculates electric potential in tissue in case of RF heating for arbitrary boundary conditions.

## 1. Theoretical basis

The RF heating is a procedure where two or more electrodes are used to apply alternating electrical current to tissue, with for example one electrode applying alternating electrical potential, while the other acts as some kind of "electrical ground"

The main equation for calculating electrical potential in the tissue is:

$$\nabla[(\sigma - i\omega\epsilon)(\nabla V)]=0$$
or when considering $\epsilon = \epsilon^/ + i\epsilon^{//}$ and that $\epsilon^/$ just gives some phase shift that does not matter when we consider RF heating with typical frequency in range of 10⁶/s:
$$\nabla[(\sigma + \omega\epsilon^{//})(\nabla V)]=0$$

Here, $\sigma$ [S/m] is electrical conductivity, i is imaginary unit, $\omega$ [1/s] is the angular frequency, $V$ is electric potential, and $\epsilon$ [F/m] is the electrical permittivity of the material.

It should be noted that the electrical permittivity, denoted here as $\epsilon$, is actually product of vacuum permittivity ${\epsilon}_0$ [F/m] with approximate value of 8.85e-12 and relative permittivity ${\epsilon}_r$ [/]. so that the following relation is valid

$$D=\epsilon E={\epsilon}_0 {\epsilon}_r E$$

## 2. Usage

The usage of the project is simple: Just run one of the example files or create a new example file and then provide the necessary parameters.

In case you want the plots saved in the Images subfolder, provide filename parameter to plot_graphs function, such as in example below:

```julia
plot_graphs(material_indices, grid_params, Qel, E_new, V_new, "Example6_Franco")
```

## 3. Literature

[1] Quasi-Static Electromagnetic Dosimetry: From Basic Principles to Examples of Applications

[2] Quasi-Static Approximation Error of Electric Field
Analysis for Transcranial Current Stimulation

[3] Numerical Study of Hyper‐Thermic Laser Lipolysis With
1,064 nm Nd:YAG Laser in Human Subjects

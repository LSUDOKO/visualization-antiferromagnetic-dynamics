# AFMVisual

A Mathematica Package for visualization of two-sublattice antiferromagnetic dynamics with arbitrary inputs, including Easy-axis anisotropy fields, hard-axis anisotropy fields, Zeeman fields, spin torques, and Dzyaloshinskii–Moriya interaction (DMI).

## Overview

AFMVisual solves the Landau-Lifshitz-Gilbert (LLG) equations in the local frame of two magnetic moments m₁ and m₂, allowing for:
- Setting up custom magnetic parameters
- Finding equilibrium positions (ground states)
- Solving and visualizing eigenmodes (resonance modes)
- Visualizing antiferromagnetic dynamics with designed driving fields (e.g., Spin-Orbit Torques)

## Installation

### Option 1: System-wide Installation
Copy the `AFMVisual.wl` file to your Mathematica applications folder:
```
$UserBaseDirectory\Applications\
```
This is typically:
```
C:\Users\[YourUsername]\AppData\Roaming\Mathematica\Applications\
```

### Option 2: Local Use
Keep `AFMVisual.wl` in your working directory and load it using the full path.

## Usage

### Loading the Package
```mathematica
(* Option 1: After system-wide installation *)
Needs["AFMVisual`"]

(* Option 2: From the current directory *)
SetDirectory["PathToYourDirectory"];
Get["AFMVisual`"]

(* Alternative syntax *)
<<"AFMVisual`"
```

### Basic Workflow

1. Set up the system's magnetic parameters:
```mathematica
(* Add one easy axis along z with amplitude=1 *)
AddEasyAxis[1, {0, 0, 1}]

(* Add one DC Zeeman field along z with amplitude=10 *)
AddBFieldDC[10, {0, 0, 1}]

(* Set AFM exchange strength to 10 *)
SetExchange[10]
```

2. Display the current system configuration:
```mathematica
DispConfig[]
```

3. Find the equilibrium position (ground state):
```mathematica
FindGS[]
```
or evolve the system to equilibrium:
```mathematica
EvolveToEq[0.1, 0.001, 10000, {1, 0, 0}, {-1, 0, 0}]
```

4. Visualize eigenmodes (resonance modes):
```mathematica
PlotEigen[]
```

5. Visualize AFM dynamics with custom driving fields:
```mathematica
ωr = N[1 + Sqrt[1*(2*10 + 1)]]*γ;
FL[t_] := 0.1*{Cos[ωr*t], Sin[ωr*t], 0};
DL[t_] := {0, 0, 0};
AFMDynamics[0.01, 0.01, 5000, FL, DL, {0, 0, 1}, {0, 0, -1}]
```

## Key Functions

- `ResetAll[]` - Reset all parameters to default values
- `SetExchange[J_]` - Set AFM exchange strength
- `AddEasyAxis[Amp_, Dir_]` - Add Easy axis with magnitude and direction
- `AddHardAxis[Amp_, Dir_]` - Add Hard axis with magnitude and direction
- `AddBFieldDC[Amp_, Dir_]` - Add DC Zeeman field with magnitude and direction
- `FindGS[]` - Find lowest energy minimum state
- `PlotEigen[]` - Plot eigenmodes for current configuration
- `AFMDynamics[]` - Visualize magnetic dynamics with custom driving fields

For a complete list of functions:
```mathematica
?AFMVisual`*
```

## Examples

See the included `Examples.nb` notebook for detailed examples demonstrating:
- Spin-flop phase visualization
- Energy minima finding
- Ground state determination
- Eigenmode visualization
- AFM resonance dynamics

## Conventions

- The angular gyromagnetic ratio is set to γ = 0.176085963023 THz*rad/T
- Effective fields are in Tesla [T] units; time scale in ps (picosecond)
- Magnetic moments are dimensionless and unitary vectors

The energy functional for the AFM system:
E[m₁,m₂] = J m₁·m₂ - Kₐ(m₁·n̂ₐ)² - Kₐ(m₂·n̂ₐ)² - Kₕ(m₁·n̂ₕ)² - Kₕ(m₂·n̂ₕ)² - H₀(m₁+m₂) + D·(m₁×m₂)

## Version History

- **Version 2.1** (2025/04/16):
  - Added support for examining eigenmodes with custom equilibrium positions
  - Improved function visibility for easier debugging

- **Version 2.0** (2025/01/28):
  - Bug fixes
  - Added support for Dzyaloshinskii–Moriya interaction (DMI)

- **Version 1.0** (2024/10/28):
  - Initial release

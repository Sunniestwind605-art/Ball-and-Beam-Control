# Ball-and-Beam Cascade Control (MATLAB)

Linear cascade PID control design and simulation for a ball-on-beam system, using pole placement to meet a 2 s settling time and 5% overshoot specification.



## Results

Simulated linear closed-loop response for $x_d = -0.25$ m, starting from $x_0 = 0$ m:

| Metric | Target | Simulated |
|---|---|---|
| Percentage overshoot | 5% | 4.84% |
| 5% settling time | 2 s | 1.44 s |
| Peak beam angle | — | 7.45° |
| Peak motor torque | — | 0.674 N·m |

These are results from the linear model integrated with `ode45`. They have not been verified on the physical apparatus.

## Approach

The controller is a two-loop cascade:

- **Inner loop (beam angle):** PD controller with a double pole at $s = -10$.
- **Outer loop (ball position):** PID controller with poles at $-1.5 \pm 1.57j$ and $-10$.
- **Setpoint prefilter:** a second-order filter $F(s) = \dfrac{K_{1i}}{K_{1d}s^2 + K_{1p}s + K_{1i}}$ that cancels the outer-loop zeros. Without it, the zeros push overshoot above the value predicted by the poles alone.

The closed-loop position response becomes $\dfrac{X(s)}{X_d(s)} = \dfrac{K_{1i}}{s^3 + K_{1d}s^2 + K_{1p}s + K_{1i}}$, which has only the designed poles.

## Repository Contents

```
.
├── BallBeam_Cascade_Workflow.m   # Gain design, simulation, metrics, animation, plots
├── report/
│   ├── BallBeam_Report.tex       # LaTeX evidence report
│   └── BallBeam_Report.pdf       # Compiled report
└── figures/                      # Root loci, step response, system diagrams
```

## Requirements

- MATLAB R2019b or later (uses `tiledlayout`)
- No additional toolboxes are required. The script uses `ode45`, which is in base MATLAB.

## Usage

Run the script from MATLAB:

```matlab
BallBeam_Cascade_Workflow               % prompts for x0 and xd, default -0.25 m
BallBeam_Cascade_Workflow(0, -0.25)     % non-interactive
BallBeam_Cascade_Workflow(0, -0.25, false)  % skip the animation
```

The script prints the gains and metrics to the Command Window, shows the beam animation, then displays the root loci and position response figures.

## Limitations

- The simulation uses the linearised model with a small-angle approximation.
- Actuator limits, sensor noise, friction, and the nonlinear dynamics are not modelled.
- Results are from simulation only.

## Author

Joseph Lobishe


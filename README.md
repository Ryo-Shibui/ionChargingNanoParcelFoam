# ionChargingNanoParcelFoam

`ionChargingNanoParcelFoam` tracks charged aerosol parcels through a prescribed
flow, electric field, and positive and negative ion fields. It adds stochastic
ion attachment to the parcel models in
[reactingChargedNanoParcelFoam](https://github.com/Ryo-Shibui/reactingChargedNanoParcelFoam).
A positive ion increments the parcel charge number by one and a negative ion
decrements it by one. The Coulomb force uses the resulting signed charge on
subsequent tracking segments.

The solver uses one-way coupling. It reads the carrier fields once, leaves the
ion concentrations unchanged when ions attach, and advances the cloud in time.
It uses the drag, Brownian, electric and other particle forces selected in the
case dictionary.

## Build with OpenCFD OpenFOAM v2406

The build requires the source and compiled custom libraries of
`reactingChargedNanoParcelFoam`. The commands below use the v2406 environment
for both builds. If the three dependency libraries are already built under
`$FOAM_USER_LIBBIN`, set `RC_NANO_SOURCE` to the installed dependency solver
directory and build only this solver.

```sh
source /path/to/OpenFOAM-v2406/etc/bashrc

git clone https://github.com/Ryo-Shibui/reactingChargedNanoParcelFoam.git
(cd reactingChargedNanoParcelFoam && ./Allwmake)
export RC_NANO_SOURCE="$(cd reactingChargedNanoParcelFoam/solver && pwd)"

git clone https://github.com/Ryo-Shibui/ionChargingNanoParcelFoam.git
cd ionChargingNanoParcelFoam
./Allwmake
ionChargingNanoParcelFoam -help
```

`RC_NANO_SOURCE` must contain the dependency's `lagrangian/` directory.
If the installed solver directories are siblings and the dependency has
`lagrangian/` directly inside its solver directory, `./Allwmake` discovers it
automatically. The required libraries are `librcNanoDistribution`,
`librcNanoIntermediate`, and `librcNanoTurbulence`. The executable is installed
in `$FOAM_USER_APPBIN`.

## Run the tutorial

From the repository root, copy the tutorial to a new case directory and run:

```sh
cp -a tutorials/positiveIonAttachment /path/to/new-case
blockMesh -case /path/to/new-case
ionChargingNanoParcelFoam -case /path/to/new-case
```

The [tutorial](tutorials/positiveIonAttachment) uses a synthetic, deliberately
large positive-ion concentration so a single initially neutral 10 nm parcel
reliably acquires positive charge in a short run. This is a smoke test with
prescribed fields. Inspect the signed charge in
`0.0005/lagrangian/reactingCloud1/chargeNumber` and the coefficient table under
`postProcessing/ionCharging/reactingCloud1/`.

## Set up a case

Start with a case configured for the dependency's `reactingCloud1`. Supply the
mesh, thermodynamic fields and dictionaries, including `p`, `T`, gas
composition and `constant/g`, as required by that cloud. Set
`application ionChargingNanoParcelFoam;` in `system/controlDict`.
`carrierFieldTime "0";` is the default. Change it to select a saved directory
on the same mesh containing all four prescribed fields.

| Field | Meaning | Units |
|---|---|---|
| `U` | Carrier-gas velocity | m/s |
| `E` | Electric field vector | V/m |
| `nP` | Positive ion number density | m⁻³ |
| `nN` | Negative ion number density | m⁻³ |

For concentrations supplied in cm⁻³, multiply by `1e6` to obtain m⁻³.
For a single ion polarity, provide the other field with zero concentration.
Set the ion mobilities as positive magnitudes and specify the attachment-model
temperature in `constant/physicalProperties`. These are example values.

```foam
T   [0 0 0 1 0 0 0]   298.15;
muP [-1 0 2 0 0 1 0]  1.79e-4;
muN [-1 0 2 0 0 1 0]  2.08e-4;
```

In `constant/reactingCloud1Properties`, set `transient true;` and
`coupled false;` in `solution`. Add the following entries to the existing
`constantProperties`. The example density and permittivity are for dry NaCl.

```foam
rho0            2165;
epsilonR        5.9;
qmin            -10;
qmax            16;
diameterChange  false;
```

`epsilonR`, `qmin`, and `qmax` are required. `diameterChange` defaults to
`false`. Keep the phase-change, devolatilisation, and surface-reaction models
set to `none` for fixed-diameter tracking. Changing particle diameter is outside
the supported workflow because attachment tables are built for fixed diameters.

Use `manualInjection` or `patchInjection`, with one `fixedValue` size
distribution per diameter. Set `parcelCharge` to the signed initial integer
charge number within `[qmin,qmax]`. The default is zero. To enable electric
migration after a neutral particle attaches its first ion, include
`CoulombSep { E E; }` under `subModels/particleForces` from the start.
The electric force is zero whenever the particle is neutral. Add the
`BrownianMotion` particle force if Brownian diffusion is needed.

## Attachment and outputs

At startup, the solver computes positive and negative attachment coefficients
for every configured diameter and integer charge in `[qmin,qmax]`. The local
attachment rates are `nP*betaP` and `nN*betaN`. The solver samples a competing
exponential waiting time within each parcel tracking segment and updates the
rates after every event. It can apply multiple events in one segment. The
carrier fields and cell ion concentrations remain fixed. Check time-step
sensitivity when charge changes and electric migration are strongly coupled.

The parcel charge bounds block transitions past the table ends. Enlarge the
bounds and rerun if the endpoints are appreciably populated. Read the
coefficient values, model names, and validity flags from
`postProcessing/ionCharging/reactingCloud1/Collision_Coefficients.csv`. In a
parallel run the master writes this single case-level table. Each saved cloud
also contains the parcel `chargeNumber` field. The
[coefficient implementation notes](Reference/EQ437_UPDATE.md) describe the
attractive and non-attractive branches.

## License

GPL-3.0-or-later. See [COPYING](COPYING).

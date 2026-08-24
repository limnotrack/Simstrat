# How to build Simstrat-AED2 with FoBiS.py
In case the AED2 library is not compiled go to ´lib/libaed2´ and run:

If you also want the newer, opt-in libaed-api based AED coupling (`ModelConfig.CoupleAED`), likewise compile `lib/libaed-water`, `lib/libaed-benthic`, `lib/libaed-demo` and `lib/libaed-api`, in that order (libaed-water must be built first). In `lib/libaed-api`, also run `make obj/aed_external.o` (or `mingw32-make obj/aed_external.o`) after the default build. This is optional: Simstrat builds and runs fine with only libaed2, using the original AED2 coupling.

### Optional: full libaed-dev path (private repos, not submodules)

`libaed-riparian`, `libaed-lighting`, `libaed-dev` and `phreeqcrm` add extra AED modules (e.g. `aed_alum`, `aed_campbell`, `aed_soilbgc`, `aed_phreeqcrm`) but come from **private** AquaticEcoDynamics repositories, so unlike `libaed2`/`libaed-water`/etc. they are never vendored, submoduled, or otherwise tracked in this repo. A developer with access to those repos can opt in locally:

1. Clone `libaed-riparian`, `libaed-lighting` (as `libaed-lighting`), `libaed-dev` and `phreeqcrm` as siblings under `lib/` (they're gitignored, so this is safe — see `.gitignore`).
2. Build them in order: `libaed-riparian`, `libaed-lighting`, then `phreeqcrm` (via CMake: `cmake .. -G "MinGW Makefiles" -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release && cmake --build .`), then `libaed-dev` with `make PHREEQDIR=../phreeqcrm`.
   - `libaed-dev`'s `aed_plots.F90` needs `libgd` (`gd.h`) for live plotting, which isn't installable on a bare MinGW toolchain (no pacman/MSYS2). It's excluded from the local `libaed-dev/Makefile` via `DEFINES+=-D_NO_PLOTS_` and dropping `aed_plots.o` from `OBJS` — a guard `aed_dev.F90` already supports upstream. Nothing else in libaed-dev is modified. This only disables GLM's interactive plotting window, which Simstrat never uses.
3. Rebuild `libaed-api`'s `obj/aed_external.o` against the full set: `make AEDRIPDIR=../libaed-riparian AEDLGTDIR=../libaed-lighting AEDDEVDIR=../libaed-dev obj/aed_external.o` (also pass `AEDBENDIR=../libaed-benthic AEDDMODIR=../libaed-demo` if rebuilding the whole `libaed-api` target, not just `aed_external.o`).
4. Build Simstrat with the `release-gnu-aed-full` FoBiS mode instead of `release-gnu`: `FoBiS.py build -mode release-gnu-aed-full`. This mode additionally links `libaed-riparian`, `libaed-lighting`, `libaed-dev` and `libPhreeqcRM`, plus `-fopenmp -lstdc++` for phreeqcrm's C++/OpenMP dependencies.

If you don't have access to these private repos, ignore this section entirely — the default `release-gnu` mode builds and runs fine without them, using only the modules in `libaed-water`/`libaed-benthic`/`libaed-demo`.

~~~bash
export F90=gfortran
mingw32-make
~~~

under windows (note that the export command only works in a bash shell) or

~~~bash
export F90=gfortran
make
~~~

under linux.

Then, from this folder (`build`), run:

~~~bash
FoBiS.py -h
~~~

to get help information about its usage.
To compile Simstrat with default configuration, run:

~~~bash
FoBiS.py build
~~~

## Clean the project
To clean the main project, simply run

~~~bash
FoBiS.py clean
~~~

whilst if you need to clean also all the libraries (e.g. if you change the final target from Linux to Win), you need to run

~~~bash
FoBiS.py rule -ex purge
~~~

## Generate the code documentation
To generate the code documentation, run

~~~bash
FoBiS.py rule -ex makedoc
~~~

The generated code documentation is saved in `doc/developer/ford/ford_doc/index.htlm`


**N.B.** you can use one-line command to call the build procedure (and others) from any folder, e.g. from `tests`
~~~bash
cd ../build; FoBiS.py build; cd -
~~~

# Computational-Methods-for-Hot-Nuclear-Matter
A stream of the course "SCI 2715 – Science Research Initiative Undergraduate Research" at University of Utah.
In this git, we collect tasks, scripts, small data files, theory, and analysis related to the semester's group project.

## The project
This project is part of a course meant to give an undergraduate some understanding of
[quantum chromodynamics](https://en.wikipedia.org/wiki/Quantum_chromodynamics) (QCD), the field theory that describes
the strong nuclear force, which binds together protons and neutrons. Our tool for study is
[lattice field theory](https://en.wikipedia.org/wiki/Lattice_field_theory).

QCD is a theory describing the interactions of quarks and gluons. If we imagine the quarks to have infinite mass
(or, formally equivalently, if we imagine there to be no quarks at all)
we have a so-called "pure SU(3) theory". In this theory, when the temperature is high enough, protons and
neutrons will "melt" into their constituents, the quarks and gluons. The temperature at which this occurs
is the deconfinement temperature, $T_d$.

In this projects students will try to extract $T_d$. Even in this simplified pure SU(3) theory, this is a moderately
challenging computational endeavor, which requires cooperation among many aspiring scientists, each managing their
own modest tasks. These tasks include, but are not limited to, the following:
1. Plan which configurations shall be generated
2. Generate configurations using [SIMULATeQCD](https://github.com/LatticeQCD/SIMULATeQCD)
3. Run measurement code from SIMULATeQCD on those configurations
4. Do statistical analysis using the [AnalysisToolbox](https://github.com/LatticeQCD/AnalysisToolbox)
5. Construct new code that automatically selects $T_d$
6. Perform a continuum-limit extrapolation
7. Collect our results and interpretations in a cogent way

In the [wiki page](https://github.com/clarkedavida/Computational-Methods-for-Hot-Nuclear-Matter/wiki), we will
collect our results and give our progress. Students will report in their own words what they have learned about
certain aspects of the theory. We will record which person carried out which task.

In the [projects space](https://github.com/clarkedavida/Computational-Methods-for-Hot-Nuclear-Matter/projects?query=is%3Aopen),
we will list the tasks that need to be done, claim them, and mark them as they are completed.

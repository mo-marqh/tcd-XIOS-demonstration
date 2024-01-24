# tcd-XIOS-demonstration
Demonstration code for XML I/O server XIOS usage.

Demonstrations of using XIOS are provided with Continuous Integration testing with respect to XIOS2 trunk.

## Environments

Environments are managed, with a little complication, to enable running on scientific desktop and on Github  Continuous Integration.

There is a helper script that one can source, in the root directory, to ensure that the scientific desktop LFRic environment is loaded, using the environment variable setup for the Makefiles, which is also compatible with the Github CI runner.

`. desktopEnv`

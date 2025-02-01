# Minimal example of path tracing in C#

This is an implementation of a modified version of [smallpt](https://www.kevinbeason.com/smallpt/) in C#, which is aimed as part of the course project of UCSB 2025 winter CS 263. The equivalent F# version can be found [here](https://github.com/LeonKang130/smallpt-fsharp).

This implementation is parallelized using `System.Threading.Tasks`. For more detail, please refer to the source code in `Program.cs`.

Here is a sample output of the program with 64 spp:
![Sample - 64spp](./sample.png)
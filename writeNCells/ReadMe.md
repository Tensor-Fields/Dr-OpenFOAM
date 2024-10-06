# Install
`git clone https://github.com/Tensor-Fields/Dr-OpenFOAM.git`

# Run
```sh
cd writeNCells
wmake
```
Now, go to your OpenFOAM case (e.g. cavity) which has a mesh (e.g. by running
`blockMesh`) and run the tool
```sh
writeNCells
```

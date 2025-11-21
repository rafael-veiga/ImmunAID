# ImmunAID
Code used in analyse for the paper:

**Patients with inflammation of unknown origin phenocopy the immune presentation of adult onset Still disease**
  
## Dependences
- **Docker** (https://www.docker.com/)
- **Nextflow** (https://www.nextflow.io/)
- **Bash** (available by default on Linux and macOS; on Windows use WSL)
This pipeline uses **containerized environments** to ensure full reproducibility across different systems. By running all steps inside Docker containers, we avoid issues related to package versions and system configurations.
  
## Instruction To Run
* **1 - create containers** :
 ```bash
# bash
chmos +x ./build_containers.sh
./build_containers.sh
```
* **2 - Run Nextflow pipeline** :
```bash
nextflow run main.nf
```

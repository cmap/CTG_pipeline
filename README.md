# CTG_pipeline
Data processing pipeline for PRISM CTG screens.

To process raw data from a screen use `pipeline.sh` or `process_CTG.R`. Or using Docker run
``` bash
docker pull cmap/ctg:latest
docker run -it
  -v ~/ctg_data:/data
  -v ~/ctg_results:/resuts \
  cmap/ctg:latest \
  -f data/mergedfile.csv \
  -m data/mapping.csv \
  -r data/raw_data.csv \
  -o results \
  -p <PROJECT>
```
replacing `<PROJECT>` with the name of the project and all other paths with the path to the appropriate data.

To generate a report for a screen use `CTG_report.Rmd` changing the data path parameter to point to processed CTG data.

The Docker image is on [Docker Hub](https://hub.docker.com/repository/docker/cmap/ctg).

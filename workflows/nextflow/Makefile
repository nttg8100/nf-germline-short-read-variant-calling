.PHONY: test-e2e clean
${HOME}/.pixi/bin/pixi:
	curl -sSL https://pixi.sh/install.sh | sh

test-e2e: ${HOME}/.pixi/bin/pixi
	${HOME}/.pixi/bin/pixi run nextflow run main.nf -profile singularity,test 

clean:
	rm -rf data/GSE110004 
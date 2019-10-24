PYTHONPATH=/usr/local/lib
export PYTHONPATH

APP=python3
LOG=flowRate.txt

all: $(LOG)
%.txt: %.py outlet.dat
	$(APP) $^ > $@
outlet.dat:
	cp ../postProcessing/averageQout/0/surfaceFieldValue.dat tmp.dat
	tr -d '()' < tmp.dat > $@
clean:
	@rm -f msm* *.png *.txt *.dat
.PHONY: clean all

CC = mpicc
PROGS = mandel mandelV1

all: $(PROGS)

plot: mandel.txt
	python view.py mandel.txt

mandel.txt: 
	./mandel > mandel.txt

plotV2: mandelV1.txt
	python view.py mandelV1.txt

mandelV1.txt: 
	./mandelV1 > mandelV1.txt

clean:
	rm -f $(PROGS) mandel.txt mandelV1.txt
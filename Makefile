TARGET = ring.png with_gravity.png without_gravity.png

all : $(TARGET)

ring.png : plot.py test.csv
	python plot.py

test.csv : alpha
	./alpha > test.csv

alpha : main.cpp
	g++ main.cpp -o alpha

with_gravity.png without_gravity.png : plot1.py with_gravity.csv without_gravity.csv
	python plot1.py

with_gravity.csv without_gravity.csv : kappa
	./kappa

kappa : mm.cpp
	g++ mm.cpp -o kappa

clean:
	rm -f $(TARGET) alpha kappa test.csv without_gravity.csv with_gravity.csv

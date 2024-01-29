compile:
	cd src;make

all:
	make clear
	make compile

clean:
	cd src;make clean
	rm -f *~ */*~

clear:
	cd src;make clear
	rm -f *~ */*~ */*.pyc


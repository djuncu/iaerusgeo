# Installation

Currently, the way to compile from the iAERUS-GEO main dir (one level above this one), is:

```bash
make clear && make clean && make && ./make.sh
```

whereas ```make.sh``` is just

```bash
cd src
./makef2py.sh 64
```

I thought that until recently the regular makefile (with some changes I had made) was still working.. I can't explain why it is not currently.


# OSX

### Lib Config (Configuration Properties)

```
$ brew install libconfig
```

### GNU Scientific Library (Science Math)

```
$ brew install gsl
```

### Argp (Argument Processing)

```
$ brew install argp-standalone
```

### CMake (For Building GoogleTest)

```
$ brew install cmake
```

### GoogleTest (C++ Unit Testing)

```
$ cd /usr/local
$ sudo git clone https://github.com/google/googletest.git
$ chmod -R dmcnelis:admin googletest
$ cd googletest
$ mkdir install
$ cd install
$ cmake -DCMAKE_CXX_COMPILER="c++" -DCMAKE_CXX_FLAGS="-std=c++11 -stdlib=libc++" ../
$ make
$ make install
```

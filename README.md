### Usage

This project is header-only, so simply include BigNum.hpp header in your project and then you can use this library like the following:

```cpp
#include <stdexcept>
#include <iostream>

#include "BigNum.hpp"

using namespace bignum;

int main(){
	bigint_t x = 114514_bigint;
	bigint_t y = x << 334;
	std::cout << std::hex << y << std::endl;
	
	return 0;
}
```

You can check our [wiki](https://github.com/gnaggnoyil/bignumplusplus/wiki) for more usages.

### Platform Support

* Clang 3.7 under x86/x86_64 linux
* Gcc 5.3 under x86/x86_64 linux
* Clang with Microsoft CodeGen under x86/x86_64 Windows

In general any platform that provides 8-bit and 32-bit integers might work. Any compiler that supports C++14 will compile this library. Compilers that support part of C++14 features might also compile, except for MSVC compilers including vc140, which has wierd issues compling some SFINAE method templates.

### TODO

* a full manual and a full document
* clean up unused code
* clean up interface
* more test
* more platform support and portability
* make BigInt work on any data length (a.k.a make BigInt's length great again)
* a better simulation for ```basic_ostream::operator>>``` so that the BigInt behaves as if a native integer.
* ```pow```, ```log```, ```exp```, ...
* minus radix support for radix conversion
* big float
* optimization
* expression template pattern to avoid some unnecessary bignum calucations
* and many more needed work to do ...

### LICENSE

You can use this library under any of the following licenses:

* MIT license
* BSD 3-clause license
* LGPL v3 license or
* GPL v2 license.
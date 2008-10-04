#ifndef DEBUG_CC
#define DEBUG_CC

#include <iostream>
#include <stdarg.h>

#ifdef EBUG
#define DEBUG(method, args...) Debug d(method, args)
#define DMORE(args...)         d.out(args)
#define DCOND(cond, args...)   if (cond) d.out(args)
#else
#define DEBUG(args...) 
#define DMORE(args...) 
#define DCOND(args...) 
#endif

class Debug {
  static int depth;

  void print(char *fmt, va_list &ap) {
    while (*fmt) {
      switch(*fmt++) {
      case 's': std::cout << va_arg(ap, char *) << " "; break;
      case 'i': std::cout << va_arg(ap, int)    << " "; break;
      case 'f': std::cout << va_arg(ap, double) << " "; break;
      }
    }
    va_end(ap);
    std::cout << std::endl;
  }

  void indent() { for (int i = 1; i < depth; i++) std::cout << "  "; }

public:
  Debug(const std::string &method, char *fmt, ...) { 
    depth++;
    indent();
    std::cout << method; if (*fmt) std::cout << ": ";
    
    va_list ap;
    va_start(ap, fmt);
    print(fmt, ap);
  }

  void out(char *fmt, ...) { 
    indent();
    va_list ap;
    va_start(ap, fmt);
    print(fmt, ap);
  }

  static void nop() {}

  ~Debug() { --depth; }
};

int Debug::depth = 0;

#endif

CFLAGS = -g
CXXFLAGS = -g -I../.. -I../../platform/macos
LD           = $(CC)

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...
all: read_options

read_options: option.o any_module.o
	$(LD) $(LDFLAGS) -o $@ $^ -lstdc++ -lutility -L/home/olange/rosetta/rosetta_source/build/src/debug/macos/10.4/32/x86/gcc/
#rosetta_options.o readinp.o

clean:
	rm *.o
	rm g_fca



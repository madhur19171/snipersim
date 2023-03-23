TARGET=fft
include ../shared/Makefile.shared

$(TARGET): $(TARGET).o
	$(CC) $(TARGET).o -lm $(SNIPER_LDFLAGS) -o $(TARGET)

run_$(TARGET):
	../../run-sniper -v -n 4 -c gainestown_log --roi -- ./fft -p 4 -m 14

CLEAN_EXTRA=viz

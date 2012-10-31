
- Out-of-memory conditions are not handled as gracefully as they might be, as 
  failed mallocs result in an abort rather than a Python MemoryError. This 
  can only happen initially and not after some time simulating, as all memory
  is alloced at creation time.

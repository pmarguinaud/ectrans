
USE DUMP_MOD
#define dump(x) CALL ODUMP (__FILE__, __LINE__, #x); CALL DUMP (x); CALL CDUMP ()

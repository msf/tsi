/* log.c */

#include <stdlib.h>
#include "tsi.h"
#include "log.h"


log_t *new_log(char *filename) {
    log_t *new_log;
    
    new_log = (log_t *) tsi_malloc(sizeof(log_t));
    if (new_log == NULL) {
        printf_dbg("new_log: failed to allocate space for new log\n");
        return NULL;
    }
    /* open file */
    return new_log;
}

char *logbuf(log_t *l) {
    if (l) return l->buf;
    printf_dbg("logbuf: received NULL as parameter!\n");
    return NULL;
}

char *reset_logbuf(log_t *l) {
    if (l) {
        l->buf[0] = 0;
        return l->buf;
    }
    printf_dbg("reset_buffer: received NULL as parameter!\n");
    return NULL;    
}

int commit_log(log_t *l) {
    if (l) {
        /* write to file */
        return 1;
    }
    printf_dbg("log: received NULL as parameter!\n");
    return 0;
}

void delete_log(log_t *l) {
    if (l) {
        /* close file */
        tsi_free(l);
    } else {
        printf_dbg("delete_log: received NULL as parameter!\n");
    }
}
 
/* end of file log.c */

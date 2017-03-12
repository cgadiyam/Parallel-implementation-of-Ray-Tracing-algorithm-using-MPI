#ifndef __SLAVE_PROCESS_H__
#define __SLAVE_PROCESS_H__

#include "RayTrace.h"

void slaveMain( ConfigData *data );

void static_strips_horizontal(ConfigData* data);

void static_cycles_vertical(ConfigData* data);

void static_blocks(ConfigData* data);

void dynamic(ConfigData* data);

#endif

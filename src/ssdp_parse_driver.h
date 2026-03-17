#ifndef SSDP_PARSE_DRIVER_H
#define SSDP_PARSE_DRIVER_H

#include "ssdp_ast.h"

typedef void *yyscan_t;

enum {
	SSDP_PARSE_COMPAT_OFF = 0,
	SSDP_PARSE_COMPAT_LEGACY = 1
};

typedef int (*ssdp_dispatch_fn)(const ssdp_cmd *cmd, void *dispatch_user);

typedef struct ssdp_parse_ctx {
	const char *filename;
	int line;
	int compat_mode;
	int strict_mode;
	int scan_col;
	ssdp_cmd *cmd;
	ssdp_dispatch_fn dispatch;
	void *dispatch_user;
} ssdp_parse_ctx;

int ssdp_parse_line(ssdp_parse_ctx *ctx, const char *line);
int ssdp_parse_and_dispatch_line(ssdp_parse_ctx *ctx, const char *line);

#endif

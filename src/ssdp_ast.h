#ifndef SSDP_AST_H
#define SSDP_AST_H

#include <stddef.h>

typedef struct ssdp_arg {
	char *key;
	char *value;
	int has_key;
	int line;
	int col;
} ssdp_arg;

typedef struct ssdp_cmd {
	char *name;
	ssdp_arg *args;
	size_t argc;
	size_t cap;
	int line;
	int col;
} ssdp_cmd;

ssdp_arg *ssdp_arg_new_take(char *key, char *value, int has_key, int line, int col);
void ssdp_arg_free(ssdp_arg *arg);

ssdp_cmd *ssdp_cmd_new_take(char *name, int line, int col);
int ssdp_cmd_add_arg_take(ssdp_cmd *cmd, ssdp_arg *arg);
void ssdp_cmd_free(ssdp_cmd *cmd);

int ssdp_cmd_get(const ssdp_cmd *cmd, const char *key, const char **value);
int ssdp_cmd_getn(const ssdp_cmd *cmd, const char *prefix, int idx, const char **value);
int ssdp_cmd_has_duplicate_key(const ssdp_cmd *cmd, const char *key);
char *ssdp_cmd_render_legacy_args(const ssdp_cmd *cmd);

char *ssdp_unquote_token(const char *token);

#endif

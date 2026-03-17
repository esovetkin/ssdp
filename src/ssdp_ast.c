#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ssdp_ast.h"

ssdp_arg *ssdp_arg_new_take(char *key, char *value, int has_key, int line, int col)
{
	ssdp_arg *arg = malloc(sizeof(*arg));
	if (NULL == arg) {
		free(key);
		free(value);
		return NULL;
	}

	arg->key = key;
	arg->value = value;
	arg->has_key = has_key;
	arg->line = line;
	arg->col = col;
	return arg;
}

void ssdp_arg_free(ssdp_arg *arg)
{
	if (NULL == arg)
		return;
	free(arg->key);
	free(arg->value);
	free(arg);
}

ssdp_cmd *ssdp_cmd_new_take(char *name, int line, int col)
{
	ssdp_cmd *cmd = malloc(sizeof(*cmd));
	if (NULL == cmd) {
		free(name);
		return NULL;
	}

	cmd->name = name;
	cmd->args = NULL;
	cmd->argc = 0;
	cmd->cap = 0;
	cmd->line = line;
	cmd->col = col;
	return cmd;
}

int ssdp_cmd_add_arg_take(ssdp_cmd *cmd, ssdp_arg *arg)
{
	ssdp_arg *tmp;

	if (NULL == cmd || NULL == arg)
		return -1;

	if (cmd->argc == cmd->cap) {
		size_t ncap = cmd->cap ? (2 * cmd->cap) : 8;
		tmp = realloc(cmd->args, ncap * sizeof(*tmp));
		if (NULL == tmp)
			return -1;
		cmd->args = tmp;
		cmd->cap = ncap;
	}

	cmd->args[cmd->argc] = *arg;
	cmd->argc++;
	free(arg);
	return 0;
}

void ssdp_cmd_free(ssdp_cmd *cmd)
{
	size_t i;

	if (NULL == cmd)
		return;

	for (i = 0; i < cmd->argc; ++i) {
		free(cmd->args[i].key);
		free(cmd->args[i].value);
	}
	free(cmd->args);
	free(cmd->name);
	free(cmd);
}

int ssdp_cmd_get(const ssdp_cmd *cmd, const char *key, const char **value)
{
	size_t i;

	if (NULL == cmd || NULL == key || NULL == value)
		return 0;

	for (i = 0; i < cmd->argc; ++i) {
		if (!cmd->args[i].has_key)
			continue;
		if (0 == strcmp(cmd->args[i].key, key)) {
			*value = cmd->args[i].value;
			return 1;
		}
	}
	return 0;
}

int ssdp_cmd_getn(const ssdp_cmd *cmd, const char *prefix, int idx, const char **value)
{
	char *key;
	int n, ret;

	if (NULL == prefix || NULL == value)
		return 0;

	n = snprintf(NULL, 0, "%s%d", prefix, idx);
	if (n < 0)
		return 0;

	key = malloc((size_t)n + 1);
	if (NULL == key)
		return 0;

	snprintf(key, (size_t)n + 1, "%s%d", prefix, idx);
	ret = ssdp_cmd_get(cmd, key, value);
	free(key);
	return ret;
}

int ssdp_cmd_has_duplicate_key(const ssdp_cmd *cmd, const char *key)
{
	size_t i, found = 0;

	if (NULL == cmd || NULL == key)
		return 0;

	for (i = 0; i < cmd->argc; ++i) {
		if (!cmd->args[i].has_key)
			continue;
		if (0 == strcmp(cmd->args[i].key, key)) {
			++found;
			if (found > 1)
				return 1;
		}
	}

	return 0;
}

char *ssdp_cmd_render_legacy_args(const ssdp_cmd *cmd)
{
	size_t i, len = 0, at = 0;
	char *out;

	if (NULL == cmd || 0 == cmd->argc) {
		out = malloc(1);
		if (out)
			out[0] = '\0';
		return out;
	}

	for (i = 0; i < cmd->argc; ++i) {
		if (cmd->args[i].has_key)
			len += strlen(cmd->args[i].key) + 1;
		len += strlen(cmd->args[i].value);
		if (i + 1 < cmd->argc)
			len += 1;
	}

	out = malloc(len + 1);
	if (NULL == out)
		return NULL;

	for (i = 0; i < cmd->argc; ++i) {
		if (cmd->args[i].has_key) {
			size_t lk = strlen(cmd->args[i].key);
			memcpy(out + at, cmd->args[i].key, lk);
			at += lk;
			out[at++] = '=';
		}

		{
			size_t lv = strlen(cmd->args[i].value);
			memcpy(out + at, cmd->args[i].value, lv);
			at += lv;
		}

		if (i + 1 < cmd->argc)
			out[at++] = ' ';
	}
	out[at] = '\0';
	return out;
}

char *ssdp_unquote_token(const char *token)
{
	size_t i, j, n;
	char q;
	char *out;

	if (NULL == token)
		return NULL;

	n = strlen(token);
	if (n < 2)
		return strdup(token);

	q = token[0];
	if (!((q == '\'' || q == '"') && token[n - 1] == q))
		return strdup(token);

	out = malloc(n - 1);
	if (NULL == out)
		return NULL;

	j = 0;
	for (i = 1; i < n - 1; ++i) {
		if (token[i] == '\\' && i + 1 < n - 1) {
			++i;
			out[j++] = token[i];
		} else {
			out[j++] = token[i];
		}
	}
	out[j] = '\0';
	return out;
}

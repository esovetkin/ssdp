#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "ssdp_parse_driver.h"
#include "ssdp_parser.tab.h"

typedef struct yy_buffer_state *YY_BUFFER_STATE;

int yylex_init(yyscan_t *scanner);
int yylex_destroy(yyscan_t scanner);
YY_BUFFER_STATE yy_scan_string(const char *str, yyscan_t scanner);
void yy_delete_buffer(YY_BUFFER_STATE b, yyscan_t scanner);

static int is_blank_line(const char *line)
{
	while (*line) {
		if (!isspace((unsigned char)*line))
			return 0;
		++line;
	}
	return 1;
}

int ssdp_parse_line(ssdp_parse_ctx *ctx, const char *line)
{
	yyscan_t scanner;
	YY_BUFFER_STATE buf;
	int ret;

	if (NULL == ctx || NULL == line)
		return -1;

	if (line[0] == '#') {
		ctx->cmd = NULL;
		return 0;
	}

	if (is_blank_line(line)) {
		ctx->cmd = NULL;
		return 0;
	}

	ctx->cmd = NULL;
	ctx->scan_col = 1;
	if (yylex_init(&scanner))
		return -1;

	buf = yy_scan_string(line, scanner);
	if (NULL == buf) {
		yylex_destroy(scanner);
		return -1;
	}

	ret = yyparse(scanner, ctx);

	yy_delete_buffer(buf, scanner);
	yylex_destroy(scanner);

	if (ret)
		return -1;

	if (NULL == ctx->cmd)
		return 0;

	return 1;
}

int ssdp_parse_and_dispatch_line(ssdp_parse_ctx *ctx, const char *line)
{
	int ret;

	ret = ssdp_parse_line(ctx, line);
	if (ret <= 0)
		return ret;

	if (NULL == ctx->dispatch)
		return 0;

	ret = ctx->dispatch(ctx->cmd, ctx->dispatch_user);
	ssdp_cmd_free(ctx->cmd);
	ctx->cmd = NULL;
	return ret;
}

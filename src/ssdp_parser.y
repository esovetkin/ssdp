%{
#include <stdio.h>
#include <stdlib.h>
#include "ssdp_parse_driver.h"
#include "ssdp_ast.h"
%}

%define api.pure full
%define parse.error detailed
%locations

%parse-param { yyscan_t scanner } { ssdp_parse_ctx *ctx }
%lex-param   { yyscan_t scanner } { ssdp_parse_ctx *ctx }

%union {
	char *str;
	ssdp_arg *arg;
	ssdp_cmd *cmd;
}

%code provides {
	int yylex(YYSTYPE *yylval_param, YYLTYPE *yylloc_param, yyscan_t scanner, ssdp_parse_ctx *ctx);
	void yyerror(YYLTYPE *loc, yyscan_t scanner, ssdp_parse_ctx *ctx, const char *msg);
}

%token <str> T_WORD
%type  <arg> item
%type  <cmd> command

%destructor { free($$); } <str>
%destructor { ssdp_arg_free($$); } <arg>
%destructor { ssdp_cmd_free($$); } <cmd>

%%
input:
	/* empty */ { ctx->cmd = NULL; }
	| command { ctx->cmd = $1; }
	;

command:
	T_WORD {
		$$ = ssdp_cmd_new_take($1, @1.first_line, @1.first_column);
		if (NULL == $$)
			YYERROR;
	}
	| command item {
		if (ssdp_cmd_add_arg_take($1, $2)) {
			ssdp_arg_free($2);
			ssdp_cmd_free($1);
			YYERROR;
		}
		$$ = $1;
	}
	;

item:
	T_WORD '=' T_WORD {
		$$ = ssdp_arg_new_take($1, $3, 1, @1.first_line, @1.first_column);
		if (NULL == $$) {
			free($1);
			free($3);
			YYERROR;
		}
	}
	| T_WORD {
		$$ = ssdp_arg_new_take(NULL, $1, 0, @1.first_line, @1.first_column);
		if (NULL == $$) {
			free($1);
			YYERROR;
		}
	}
	;
%%

void yyerror(YYLTYPE *loc, yyscan_t scanner, ssdp_parse_ctx *ctx, const char *msg)
{
	int line = 1;
	int col = 1;
	const char *fn = "<stdin>";
	(void) scanner;

	if (ctx) {
		if (ctx->line > 0)
			line = ctx->line;
		if (ctx->filename)
			fn = ctx->filename;
	}

	if (loc) {
		if (loc->first_column > 0)
			col = loc->first_column;
	}

	fprintf(stderr, "Parse error at %s:%d:%d: %s\n", fn, line, col, msg);
}

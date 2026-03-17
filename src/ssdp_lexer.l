%option reentrant noyywrap nounput noinput bison-bridge bison-locations

%{
#include <stdlib.h>
#include "ssdp_parse_driver.h"
#include "ssdp_parser.tab.h"
#include "ssdp_ast.h"

#define YY_DECL int yylex(YYSTYPE *yylval_param, YYLTYPE *yylloc_param, yyscan_t yyscanner, ssdp_parse_ctx *ctx)

#define YY_USER_ACTION \
	yylloc_param->first_line = 1; \
	yylloc_param->last_line = 1; \
	yylloc_param->first_column = ctx->scan_col; \
	yylloc_param->last_column = ctx->scan_col + (int)yyleng - 1; \
	ctx->scan_col += (int)yyleng;
%}

%%
[ \t\r]+                         ;
"="                              { return '='; }
\"([^\"\\\n]|\\.)*\"             { yylval_param->str = ssdp_unquote_token(yytext); return T_WORD; }
\'([^\'\\\n]|\\.)*\'             { yylval_param->str = ssdp_unquote_token(yytext); return T_WORD; }
[^ \t\r\n=]+                     { yylval_param->str = strdup(yytext); return T_WORD; }
\n                               { return 0; }
%%

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../src/ssdp_parse_driver.h"

static int parse_cmd(const char *line, ssdp_cmd **out)
{
	ssdp_parse_ctx ctx;
	int ret;

	memset(&ctx, 0, sizeof(ctx));
	ctx.compat_mode = SSDP_PARSE_COMPAT_LEGACY;
	ret = ssdp_parse_line(&ctx, line);
	if (ret > 0)
		*out = ctx.cmd;
	else
		*out = NULL;
	return ret;
}

static void test_simple_command(void)
{
	ssdp_cmd *cmd = NULL;
	assert(1 == parse_cmd("config_sky C=C N=21", &cmd));
	assert(cmd);
	assert(0 == strcmp(cmd->name, "config_sky"));
	ssdp_cmd_free(cmd);
}

static void test_numbered_args(void)
{
	const char *v = NULL;
	ssdp_cmd *cmd = NULL;
	assert(1 == parse_cmd("write_h5 a0=A a1=B file=x.h5", &cmd));
	assert(cmd);
	assert(ssdp_cmd_getn(cmd, "a", 0, &v));
	assert(0 == strcmp(v, "A"));
	assert(ssdp_cmd_getn(cmd, "a", 1, &v));
	assert(0 == strcmp(v, "B"));
	ssdp_cmd_free(cmd);
}

static void test_quoted_value(void)
{
	const char *v = NULL;
	ssdp_cmd *cmd = NULL;
	assert(1 == parse_cmd("read_sky file=\"a b.h5\" dataset=sky", &cmd));
	assert(cmd);
	assert(ssdp_cmd_get(cmd, "file", &v));
	assert(0 == strcmp(v, "a b.h5"));
	ssdp_cmd_free(cmd);
}

static void test_single_quoted_value(void)
{
	const char *v = NULL;
	ssdp_cmd *cmd = NULL;
	assert(1 == parse_cmd("read_sky file='a b.h5' dataset=sky", &cmd));
	assert(cmd);
	assert(ssdp_cmd_get(cmd, "file", &v));
	assert(0 == strcmp(v, "a b.h5"));
	ssdp_cmd_free(cmd);
}

static void test_value_characters(void)
{
	const char *v = NULL;
	ssdp_cmd *cmd = NULL;
	assert(1 == parse_cmd("read_h5 file=/tmp/a.b/${x}:1 dataset=d", &cmd));
	assert(cmd);
	assert(ssdp_cmd_get(cmd, "file", &v));
	assert(0 == strcmp(v, "/tmp/a.b/${x}:1"));
	ssdp_cmd_free(cmd);
}

static void test_malformed_assignments(void)
{
	ssdp_cmd *cmd = NULL;
	assert(-1 == parse_cmd("config_sky x=", &cmd));
	assert(NULL == cmd);
	assert(-1 == parse_cmd("config_sky =x", &cmd));
	assert(NULL == cmd);
	assert(-1 == parse_cmd("config_sky x==1", &cmd));
	assert(NULL == cmd);
}

static void test_comment_and_blank_line(void)
{
	ssdp_cmd *cmd = NULL;
	assert(0 == parse_cmd("# comment", &cmd));
	assert(NULL == cmd);
	assert(0 == parse_cmd("    ", &cmd));
	assert(NULL == cmd);
}

static void test_duplicate_key_policy(void)
{
	const char *v = NULL;
	ssdp_cmd *cmd = NULL;
	assert(1 == parse_cmd("make_scalar x=a val=1 val=2", &cmd));
	assert(cmd);
	assert(ssdp_cmd_get(cmd, "val", &v));
	assert(0 == strcmp(v, "1"));
	assert(ssdp_cmd_has_duplicate_key(cmd, "val"));
	ssdp_cmd_free(cmd);
}

static void test_exit(void)
{
	ssdp_cmd *cmd = NULL;
	assert(1 == parse_cmd("exit", &cmd));
	assert(cmd);
	assert(0 == strcmp(cmd->name, "exit"));
	ssdp_cmd_free(cmd);
}

int main(void)
{
	test_simple_command();
	test_numbered_args();
	test_quoted_value();
	test_single_quoted_value();
	test_value_characters();
	test_malformed_assignments();
	test_comment_and_blank_line();
	test_duplicate_key_policy();
	test_exit();
	printf("ssdptest_parser: PASS\n");
	return 0;
}

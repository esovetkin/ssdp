#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#ifdef _WIN64
#include "winunistd.h"
#else
#include <unistd.h>
#endif /*_WIN64*/
#include <locale.h>
#include <readline/readline.h>
#include <readline/history.h>
#include "parser.h"



int running;
int sigwinch_received;
const char *prompt = ">> ";

/* Handle SIGWINCH and window size changes when readline is not active and
   reading a character. */
/* SIGWINCH is not present on windows, I just ignore this part then ... */
#ifndef _WIN64
static void sighandler (int sig)
{
	sigwinch_received = 1;
}
#endif /*_WIN64*/


char ** keyword_completion(const char *text, int start, int end)
{
  char **matches;
  matches = (char **)NULL;
  /* If this word is at the start of the line, then it is a command
     to complete.  Otherwise it is the name of a file in the current
     directory. */
  if (start == 0)
    matches = rl_completion_matches (text, keyword_generator);
  return (matches);
}
/* Callback function called for each line when accept-line executed, EOF
   seen, or EOF character read.  This sets a flag and returns; it could
   also call exit(3). */
static void cb_linehandler (char *line)
{
	/* Can use ^D (stty eof) or `exit' to exit. */
	if (line == NULL || strcmp (line, "exit") == 0)
    {
		if (line == 0)
			printf ("\n");
		printf ("bye\n");
      /* This function needs to be called to reset the terminal settings,
         and calling it from the line handler keeps one extra prompt from
         being displayed. */
		rl_callback_handler_remove ();

		running = 0;
	}
	else
	{
		if (*line)
			add_history (line);
		
		ParseComm(line);
		free (line);
	}
}

void shell()
{
	int r;
    rl_attempted_completion_function = keyword_completion;
	setlocale (LC_ALL, "");
#ifndef _WIN64
	signal (SIGWINCH, sighandler);
#endif /*_WIN64*/
	rl_callback_handler_install (prompt, cb_linehandler);
	running = 1;
	while (running)
    {
		r = isatty(STDIN_FILENO);
		if (r != 1)
		{
			fprintf(stderr,"Warning: Problems reading input");
			rl_callback_handler_remove ();
			break;
		}
		if (sigwinch_received)
		{
			rl_resize_terminal ();
			sigwinch_received = 0;
		}
		if (r == 1)
			rl_callback_read_char ();
	}
}

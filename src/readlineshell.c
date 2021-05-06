#include <stdlib.h>
#include <stdio.h>
#ifdef _WIN32
void shell()
{
	fprintf(stderr,"Sorry, no interactive shell on windows\n");
}
#else
#include <string.h>
#include <unistd.h>
#include <locale.h>

/* Used for select(2) */
#include <sys/types.h>
#include <sys/select.h>

#include <signal.h>

#include <stdio.h>
#include <errno.h>

/* Standard readline include files. */
#include <readline/readline.h>
#include <readline/history.h>


/* local includes */
#include "parser.h"


static void cb_linehandler (char *);
static void sighandler (int);


int running;
int sigwinch_received;
const char *prompt = ">> ";

/* Handle SIGWINCH and window size changes when readline is not active and
   reading a character. */
static void sighandler (int sig)
{
	sigwinch_received = 1;
}


char ** keyword_completion(const char *text, int start, int end)
{
    rl_attempted_completion_over = 1;
    return rl_completion_matches(text, keyword_generator);
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
	fd_set fds;
	int r;
    rl_attempted_completion_function = keyword_completion;

	/* Set the default locale values according to environment variables. */
	setlocale (LC_ALL, "");


	/* Handle window size changes when readline is not active and reading
		characters. */
	signal (SIGWINCH, sighandler);

	  /* Install the line handler. */
	rl_callback_handler_install (prompt, cb_linehandler);

	  /* Enter a simple event loop.  This waits until something is available
     to read on readline's input stream (defaults to standard input) and
     calls the builtin character read callback to read it.  It does not
     have to modify the user's terminal settings. */
	running = 1;
	while (running)
    {
		FD_ZERO (&fds);
		FD_SET (fileno (rl_instream), &fds);    

		r = select (FD_SETSIZE, &fds, NULL, NULL, NULL);
		if (r < 0 && errno != EINTR)
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
		if (r < 0)
			continue;     

		if (FD_ISSET (fileno (rl_instream), &fds))
			rl_callback_read_char ();
	}
}
#endif /*_WIN32 */

#include<string>
class align_wrap
{
	public:
	/* 0 = ran to completion, 1 = bad arguments (usage was printed). */
	int status;
	~align_wrap();
	align_wrap(char *argv[]);
};

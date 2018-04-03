#include <stdio.h>
#include <string.h>

int main ( int argc, char *argv[] )
{
	char filename[256];
	char filename_out[256];
    	char buf[2000];
	if (argc==2)
	{
		strcpy(filename,argv[1]);
		strcpy(filename_out,argv[1]);
		strcat(filename_out,".cmp1");
		FILE *fp = fopen (filename, "r");
		FILE *fop = fopen (filename_out, "w");
		if(fp!=NULL)
		{
			char * line = NULL;
			size_t len = 0;
			ssize_t read;
			while ((read = getline(&line, &len, fp)) != -1) 
			{
				if (line[0]!='>')
				{
					int i=0;
					for(i=strlen(line);i>0 && (line[i]=='\n' || line[i]=='\r' || line[i]==' ' || line[i]=='\t' || line[i]=='\0');i--) { line[i]='\0'; }
					fprintf(fop,"%s", line);
				}
			}
		}
		fclose (fp);
	}
    return 0;
}

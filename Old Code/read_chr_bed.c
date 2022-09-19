#include <stdio.h>
#include <string.h>

int main ( int argc, char *argv[] )
{
	char filename[512];
	char filename_out[512];
    char chr[256];
    char chr_dir[2000];
    char name[2000];
    char qual[256];
    char strand[256];
    char lengths[200000];
    char offsets[200000];
    char buf[500000];
	if (argc>=2)
	{
		strcpy(filename,argv[1]);
		strcpy(filename_out,argv[1]);
		if (argc>=3) { strcpy(chr_dir,argv[2]); } else { strcpy(chr_dir,"/ifs/data/proteomics/tcga/databases/genome_human"); } 
		strcat(filename_out,".dna");
		FILE *fp = fopen (filename, "r");
		FILE *fop = fopen (filename_out, "w");
		char * line = NULL;
		size_t len = 0;
		ssize_t read;
		while ((read = getline(&line, &len, fp)) != -1) 
		{
			int i=0;
			int j=0;
			for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { chr[j]=line[i]; }
			chr[j]='\0';
			if (line[i]!='\0')
			{
				i++;
				int start=atoi(line+i);
				for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
				if (line[i]!='\0')
				{
					i++;
					for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
					if (line[i]!='\0')
					{
						i++;
						for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { name[j]=line[i]; }
						name[j]='\0';
						if (line[i]!='\0')
						{
							i++;
							for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { qual[j]=line[i]; }
							qual[j]='\0';
							if (line[i]!='\0')
							{
								i++;
								for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { strand[j]=line[i]; }
								strand[j]='\0';
								if (line[i]!='\0')
								{
									i++;
									for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
									if (line[i]!='\0')
									{
										i++;
										for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
										if (line[i]!='\0')
										{
											i++;
											for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
											if (line[i]!='\0')
											{
												i++;
												for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
												if (line[i]!='\0')
												{
													i++;
													for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { lengths[j]=line[i]; }
													lengths[j]='\0';
													if (line[i]!='\0')
													{
														i++;
														for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { offsets[j]=line[i]; }
														offsets[j]='\0';
														if (line[i]!='\0')
														{
															i++;
															fprintf(fop,"%s",line);
															fprintf(fop,">%s (MAP:%s:%d%s %s %s)\n",name,chr,start,strand,lengths,offsets);
															int k=0;
															int l=0;
															int m=0;
															char filename_chr[256];
															int padding=600;
															strcpy(filename_chr,chr_dir);
															strcat(filename_chr,"/");
															strcat(filename_chr,chr);
															strcat(filename_chr,".fa.cmp1");
															FILE *fp_chr = fopen(filename_chr, "r");
															if (fp_chr!=NULL)
															{
																buf[0]='\0';
																fseek (fp_chr, start-padding, SEEK_SET);
																fgets(buf, padding+1, fp_chr);
																if (buf!=NULL)
																{
																	fprintf(fop,"%s\t-1\t%s\t%d\t%d\t%s\t0\t%s\n",name,chr,start-padding,padding,strand,buf);
																}
																int start_=0;
																int length=0;
																for(k=0,l=0,m=0;offsets[l]!='\0';k++)
																{
																	if (isdigit(offsets[l]))
																	{
																		start_=start+atoi(offsets+l);
																		if (isdigit(lengths[m]))
																		{
																			length=atoi(lengths+m);
																			buf[0]='\0';
																			fseek (fp_chr, start_, SEEK_SET);
																			if (length<500000)
																			{
																				fgets(buf, length+1, fp_chr);
																				if (buf!=NULL)
																				{
																					fprintf(fop,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\n",name,k,chr,start_,length,strand,qual,buf);
																				}
																			} else { fprintf(fop,"%s\t%d\t%s\t%d\t%d\t%s\t%s\tError Reading Sequence\n",name,k,chr,start_,length,strand,qual); }
																		}
																	}
																	for(;offsets[l]!=',' && offsets[l]!='\0';l++) { ; }
																	if (offsets[l]==',') { l++; }
																	for(;lengths[m]!=',' && lengths[m]!='\0';m++) { ; }
																	if (lengths[m]==',') { m++; }
																}
																buf[0]='\0';
																fseek (fp_chr, start_+length, SEEK_SET);
																fgets(buf, padding+1, fp_chr);
																if (buf!=NULL)
																{
																	fprintf(fop,"%s\t+1\t%s\t%d\t%d\t%s\t0\t%s\n",name,chr,start_+length,padding,strand,buf);
																}
																fclose(fp_chr);
															} else { printf("Error opening: %s (sequence likely not on a 'conventional' chromosome)\n\n",filename_chr); }
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
    return 0;
}

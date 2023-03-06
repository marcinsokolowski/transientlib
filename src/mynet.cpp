#include <arpa/inet.h>
#include <netinet/in.h>
#include <unistd.h>
#include <netdb.h>
#include <errno.h>

#include "mynet.h"



int get_my_ip( mystring& szIP )
{
	/* have you ever seen a hostname longer than a screen (80cols)?*/
	char name[512]; /*store my hostname*/
	struct hostent * hostent_ptr;
	int ret;
	ret = gethostname (name, 80);
	if(ret == -1) {	
		printf ("ERROR gethostname() failed, Errno=%d \nDescription: %s\n", errno,strerror(errno));
		return 0;
	}


	hostent_ptr = gethostbyname(name);
	if(hostent_ptr == NULL)	
	{
		printf ("ERROR gethostbyname() failed, h_errno=%d \nDescription: %s\n",h_errno, hstrerror(h_errno));
		return 0;
	}

	/*h_addr_list contains IPs of this host in network byte order */

	in_addr_t my_ip = ((struct in_addr *)hostent_ptr->h_addr_list[0])->s_addr; /*get the first IP.*/		

	struct in_addr temp;
                                                                                
	temp.s_addr = my_ip;                                                                                
	szIP = inet_ntoa(temp);

	return 1;	
}

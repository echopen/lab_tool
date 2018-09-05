/*----
Author: Jérôme Dubois
Contributors:
Version: 1.1
Date: 08/2018
Descritption:

A C librairy create to manage server and client using socket via TCP protocol. Transmited data can be char or int16_t.

Initialy create to make a server on a RedPitaya defined as access point and receive the data on PC client. Server can handle multiple connection, and send the same data at each client.

----*/

#include<stdio.h> //for printf
#include<stdlib.h> //for exit()
#include<stdint.h> //int16_t etc
#include<sys/socket.h> //for socket
#include<sys/types.h> //utile?
#include<arpa/inet.h> //for inet_addr
#include<errno.h> //for error
#include<unistd.h> //for close()
#include<netdb.h> //for gethostbyname
#include<pthread.h> //for thread

#include"TCP_API.hpp"

pthread_t TCP_server_thread;

void init_struct_client(client* client_list, unsigned int Nmax)
{
	unsigned int i=0;
	client_list->Nmax=Nmax;
	client_list->NbClient=0;
	client_list->sock_client=(SOCKET *)malloc(Nmax*sizeof(SOCKET));
	client_list->sin_client=(SOCKADDR_IN *)malloc(Nmax*sizeof(SOCKADDR_IN));
	client_list->id_client=(unsigned int *)malloc(Nmax*sizeof(unsigned int));
	for (i=0 ; i<Nmax ; i++){client_list->id_client[i]=0;}
}

void clear_struct_client(client* client_list)
{
	free(client_list->sock_client);
	free(client_list->sin_client);
	free(client_list->id_client);
}

void add_client(client* client_list, SOCKET sock_server)
{
		unsigned int socklen=sizeof(SOCKADDR_IN);
		SOCKET tmp;
		//static unsigned int id=0;
		int Nset=6;
		char buffer[Nset];

		buffer[0]=(char)r0;
		buffer[1]=(char)rf;
		buffer[2]=(char)decimation;
		buffer[3]=(char)Nline;
		buffer[4]=(char)sector;
		buffer[5]=(char)mode_RP;

		if (client_list->NbClient<=client_list->Nmax)
		{
			tmp=accept(sock_server,(SOCKADDR *)&client_list->sin_client[client_list->NbClient],&socklen);
			client_list->sock_client[client_list->NbClient]=tmp; //strangely we have to pass by a temporary variable SOCKET or their is a problem when transmetting a socket after we clear a client
			if (client_list->sock_client[client_list->NbClient]==INVALID_SOCKET)
			{
				perror("accept()");
				exit(errno);
			}
			else
			{
				send_TCP_server(client_list, buffer, Nset, client_list->NbClient);
				client_list->NbClient+=1;
				printf("Client number %i on %i, connected on socket = %i\n",client_list->NbClient,client_list->Nmax,client_list->sock_client[client_list->NbClient-1]);
			}
		}
		//else{printf("too many clients\n");}
}

void clear_client(client* client_list, unsigned int id)
{
	unsigned int i;
	client client_temp;//{.NbClient=0, .Nmax=0, .sock_client=NULL, .sin_client=NULL, .id_client=NULL};
	init_struct_client(&client_temp,client_list->Nmax);
	if (client_list->NbClient>1)
	{
		//client_temp.NbClient=client_list->NbClient-1;
		if (id==client_list->NbClient-1)
		{
			for (i=0 ; i<id ; i++)
			{
				client_temp.sock_client[i]=client_list->sock_client[i];
				client_temp.sin_client[i]=client_list->sin_client[i];
				client_temp.id_client[i]=client_list->id_client[i];
			}
		}
		else
		{
			for (i=0 ; i<id ; i++)
			{
				client_temp.sock_client[i]=client_list->sock_client[i];
				client_temp.sin_client[i]=client_list->sin_client[i];
				client_temp.id_client[i]=client_list->id_client[i];
			}
			for (i=id+1 ; i<client_list->NbClient ; i++)
			{
				client_temp.sock_client[i-1]=client_list->sock_client[i];
				client_temp.sin_client[i-1]=client_list->sin_client[i];
				client_temp.id_client[i-1]=client_list->id_client[i];
			}
		}
	}

	for (i=0 ; i<client_list->NbClient-1 ; i++)
	{
		client_list->sock_client[i]=client_temp.sock_client[i];
		client_list->sin_client[i]=client_temp.sin_client[i];
	}
	client_list->NbClient-=1;
}

void init_TCP_client(SOCKET* sock, const char* IP, int Port)
{
	SOCKADDR_IN sclient;//={0};

	//Create socket
	(*sock)=socket(AF_INET, SOCK_STREAM, 0);
	if ((*sock)==INVALID_SOCKET)
	{
		perror("socket()");
		exit(errno);
	}

	//Socket info
	sclient.sin_addr.s_addr=inet_addr(IP);
	sclient.sin_family=AF_INET;
	sclient.sin_port=htons(Port);

	//Connection to server
	if (connect((*sock), (SOCKADDR *) &sclient, sizeof(SOCKADDR))==SOCKET_ERROR)
	{
		perror("connect()");
		exit(errno);
	}
	printf("Connected\n");
}

int int_converter(char x)
{
	int y=0;
	y=(int)x+256;
	y=y%256;
	return y;
}

void get_RP_settings(SOCKET *sock)
{
	int Nset=6;
	char tmp[Nset];
	receive_TCP_client(sock, tmp, Nset);
	r0=(float)int_converter(tmp[0]);
	rf=(float)int_converter(tmp[1]);
	decimation=int_converter(tmp[2]);
	Nline=int_converter(tmp[3]);
	sector=(double)int_converter(tmp[4]);
	mode_RP=int_converter(tmp[5]);
}

void init_TCP_server(SOCKET* sock, int Port,client* client_list, unsigned int MaxClient)
{
	SOCKADDR_IN server;//={0};

	//Create socket
	(*sock)=socket(AF_INET, SOCK_STREAM,0);
	if ((*sock)==INVALID_SOCKET)
	{
		perror("socket()");
		exit(errno);
	}

	//Socket info
	server.sin_family=AF_INET;
	server.sin_port=htons(Port);
	server.sin_addr.s_addr=INADDR_ANY; //allow all IP to access to the server

	//Link communication point
	if (bind((*sock), (SOCKADDR *)&server, sizeof(SOCKADDR))==SOCKET_ERROR)
	{
		perror("bind()");
		exit(errno);
	}

	//Connexion
	if (listen((*sock), MaxClient)==SOCKET_ERROR)
	{
		perror("listen()");
		exit(errno);
	}

	init_struct_client(client_list, MaxClient);

}

void *TCP_server_routine(void* p_data)
{
	/*server_info* serv_info=p_data;

	while (1)
	{
		add_client(serv_info->client_list,serv_info->sock);
		printf("client number %i connected\n",serv_info->client_list->NbClient);
	}*/

	return NULL;
}

void launch_server(SOCKET* sock, client* client_list)
{
	server_info* serv_info=(server_info *)malloc(sizeof(server_info));
	serv_info->sock=(*sock);
	serv_info->client_list=client_list;

	if(pthread_create(&TCP_server_thread, NULL, TCP_server_routine, serv_info)!=0)
	{
		perror("pthread_create()");
		exit(errno);
	}
}

void close_TCP_server(SOCKET* sock, client* client_list)
{
	unsigned int i=0;

	pthread_cancel(TCP_server_thread);//close thread

	for(i=0 ; i<client_list->NbClient ; i++){close(client_list->sock_client[i]);}
	close((*sock));
	clear_struct_client(client_list);
}

void close_TCP_client(SOCKET* sock)
{
	close((*sock));
}

int send_TCP_server(client* client_list, char* buffer, int buff_length, int target)
{
	unsigned int i=0;
	int err=0;

	//local variable to clear the client after sending, whereas due to reorganisation it bug if we don't clear the last client
	int tmp[client_list->Nmax+1];
	for (i=0 ; i<client_list->Nmax+1 ; i++){tmp[i]=0;}

	if (target<0)
	{
		for (i=0 ; i<client_list->NbClient ; i++)
		{
			if(send(client_list->sock_client[i], buffer, buff_length, MSG_NOSIGNAL)<0)
			{
				printf("Client %i disconnected\n",client_list->id_client[i]+1);
				//clear_client(client_list, client_list->id_client[i]);
				tmp[0]+=1;
				tmp[tmp[0]]=client_list->id_client[i];
				err=1;
			}
		}

		//Clearing disconnected clients
		for (i=tmp[0] ; i>0 ; i--)
		{
			clear_client(client_list, tmp[i]);
		}
	}
	else
	{
		if(send(client_list->sock_client[target], buffer, buff_length, MSG_NOSIGNAL)<0)
		{
			printf("Client %i disconnected\n",client_list->id_client[target]+1);
			clear_client(client_list, client_list->id_client[target]);
			err=1;
		}
	}

	return err;
}

void send_TCP_client(SOCKET* sock, char* buffer, int buff_length)
{
	if(send((*sock), buffer, buff_length, 0)<0)
	{
		perror("send()");
		exit(errno);
	}
}

int receive_TCP_server(client* client_list, char* buffer, int buff_length, int target)
{
	if (recv(client_list->sock_client[target], buffer, buff_length+1, MSG_WAITALL)==0)
	{
		printf("client disconnected\n");
		return 1;
	}
	else {return 0;}
}

int receive_TCP_client(SOCKET* sock, char* buffer, int buff_length)
{
	if (recv((*sock), buffer, buff_length, MSG_WAITALL)==0)
	{
		printf("server closed\n");
		return 1;
	}
	else {return 0;}
}

int send_int16_TCP_server(client* client_list, int16_t *buffer, int buff_length, int target)
{
	unsigned int i=0;
	int err=0;

	//local variable to clear the client after sending, whereas due to reorganisation it bug if we don't clear the last client
	int tmp[client_list->Nmax+1];
	for (i=0 ; i<client_list->Nmax+1 ; i++){tmp[i]=0;}

	if (target<0)
	{
		for (i=0 ; i<client_list->NbClient ; i++)
		{
			if(send(client_list->sock_client[i], (char *)buffer,2*buff_length, MSG_NOSIGNAL)<0)
			{
				printf("Client %i disconnected\n",client_list->id_client[i]+1);
				//clear_client(client_list, client_list->id_client[i]);
				tmp[0]+=1;
				tmp[tmp[0]]=client_list->id_client[i];
				err=1;
			}
		}

		//Clearing disconnected clients
		for (i=tmp[0] ; i>0 ; i--)
		{
			clear_client(client_list, tmp[i]);
		}
	}
	else
	{
		if(send(client_list->sock_client[target], (char *)buffer, 2*buff_length, MSG_NOSIGNAL)<0)
		{
			printf("Client %i disconnected\n",client_list->id_client[target]+1);
			clear_client(client_list, client_list->id_client[target]);
			err=1;
		}
	}

	return err;
}

void send_int16_TCP_client(SOCKET* sock, int16_t *buffer, int buff_length)
{
	if(send((*sock), (char *)buffer, 2*buff_length, 0)<0)
	{
		perror("send()");
		exit(errno);
	}
}

int receive_int16_TCP_client(SOCKET* sock, int16_t *buffer, int buff_length)
{
	//int iloop=0;
	if (recv((*sock), (char *)buffer, 2*buff_length, MSG_WAITALL)==0)
	{
		printf("server closed\n");
		return 1;
	}
	else
	{
		/*if (buff_length>1)
		{
			for (iloop=0 ; iloop<5 ; iloop++)
			{
				printf("line[%i] = %i\n", iloop, buffer[iloop]);
			}
		}*/
		return 0;
	}
}

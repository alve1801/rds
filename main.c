/*
this one is supposed to calculate what a general gradient of the two reactants would look like
tbh id like to let it run overnight and see what happens
we might wanna trim it to some more relevant parts tho
also, this varies the initial reactant concentration, which we don't really care about
hm. random initial setup, and gradient of reaction parameters?
	itll require to rewrite the reaction code a bit weirdly
...uh, nope. this actually also varies the reaction parameters.

we might wanna figure out how to let structures develop on a larger scale
like, is it just about keeping step big?
yeah, no, that doesnt change shit
multiple diffusions/reaction?
we had 30, now trying 100
yep, that seems to do the trick
we can actually adjust that over the diffusion coefficents...

toroids are cool, but uuuh we might wanna do solid walls
	at least for the gradient tests
basically, reflect it back onto itself
we can modify the mod to do that
im... not sure whether this actually works?

*/

#include <stdio.h>
#include <stdlib.h>
#define size 600
#define iterations 1000
#define step 2
#define rat 1
#define da .04
#define db .02
#define F .053
#define k .06

#define scalex 16
#define scaley 16

//int mod(int a,int b){a%=b;return a<0?a+b:a;}

int mod(int a,int b){
	//printf("%i %i\n",a,b);
	if(a<0)
		return mod(-1*a,b);
	if(a>=b)
		return mod(b-a+b-1,b);
	//putchar(10);
	return a;
}

typedef struct{float a,b;} Chem;

Chem substrate[size*size];

float r(){return(float)(rand()%100)/100;}

void init(){
	for(int x=0;x<size;x++)
	for(int y=0;y<size;y++)
	substrate[x*size+y]=(Chem){r(),r()};
	//substrate[x*size+y]=(Chem){(float)x/size,(float)y/size};
}

void diffuse(){
	Chem t[size*size];float tmp=0;
	for(int x=0;x<size;x++)
	for(int y=0;y<size;y++){

		// fuckin laplacian
		t[x*size+y].a=
			substrate[x*size+y].a * -4 +
			substrate[x*size+mod(y-1,size)].a +
			substrate[mod(x-1,size)*size+y].a +
			substrate[x*size+mod(y+1,size)].a +
			substrate[mod(x+1,size)*size+y].a
			;

		t[x*size+y].b=
			substrate[x*size+y].b * -4 +
			substrate[x*size+mod(y-1,size)].b +
			substrate[mod(x-1,size)*size+y].b +
			substrate[x*size+mod(y+1,size)].b +
			substrate[mod(x+1,size)*size+y].b
			;
	}

	for(int i=0;i<size*size;i++){
		substrate[i].a+=t[i].a*da/step;
		substrate[i].b+=t[i].b*db/step;
	}
}


void react(){
	Chem t[size*size];float tmp=0;
	for(int x=0;x<size;x++)
	for(int y=0;y<size;y++){

		// fuckin laplacian
		t[x*size+y].a=
			substrate[x*size+y].a
			;

		t[x*size+y].b=
			substrate[x*size+y].b
			;

		tmp = substrate[x*size+y].a*substrate[x*size+y].b*substrate[x*size+y].b;

		//t[x*size+y].a = - tmp + F*(1-substrate[x*size+y].a);
		//t[x*size+y].b = + tmp - (F+k)*substrate[x*size+y].b;

		t[x*size+y].a = - tmp + ((float)x/size/scalex)*(1-substrate[x*size+y].a);
		t[x*size+y].b = + tmp - ((float)x/size/scalex+(float)y/size/scaley)*substrate[x*size+y].b;

	}

	for(int i=0;i<size*size;i++){
		substrate[i].a+=t[i].a/step;
		substrate[i].b+=t[i].b/step;
	}
}

void print(){
	for(int x=0;x<size;x++){for(int y=0;y<size;y++)
	printf("<%.3f %.3f>",substrate[x*size+y].a,substrate[x*size+y].b);
	putchar(10);}putchar(10);
}

void p(){
	for(int x=0;x<size;x++){
		for(int y=0;y<size;y++)putchar(" -.:=+*%#@"[(int)(substrate[x*size+y].a*10)]);
		printf(" | ");
		for(int y=0;y<size;y++)putchar(" -.:=+*%#@"[(int)(substrate[x*size+y].b*10)]);
		putchar(10);
	}
	putchar(10);
}

void pb(){
	char s[]=" `.~-':,_\"!;^+i><?)l{[}I1](|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";

	for(int i=0;i<70;i++)putchar(s[i]);putchar(10);
	for(int x=0;x<size;x++){
		for(int y=0;y<size;y++)putchar(s[(int)(substrate[x*size+y].a*70)]);
		printf(" | ");
		for(int y=0;y<size;y++)putchar(s[(int)(substrate[x*size+y].b*70)]);
		putchar(10);
	}
	putchar(10);
}

void main(){
	init();
	FILE*f;
	char fname[16];

	f=fopen("a19.ppm","r");
	fprintf(f,"P5\n%i %i\n255\n",size,size);
	for(int i=0;i<15;i++)getc(f);
	for(int i=0;i<size*size;i++)
		substrate[i].a=(float)getc(f)/255;
	fclose(f);

	f=fopen("b19.ppm","r");
	fprintf(f,"P5\n%i %i\n255\n",size,size);
	for(int i=0;i<15;i++)getc(f);
	for(int i=0;i<size*size;i++)
		substrate[i].b=(float)getc(f)/255;
	fclose(f);

	for(int frame=20;frame<50;frame++){

		for(int i=0;i<iterations;i++){
			printf("Frame %i, Iteration %i of %i\n",frame,i,iterations);
			for(int j=0;j<rat;j++)
				diffuse();
			react();
		}

		//print();

		sprintf(fname,"a%i.ppm",frame);
		f=fopen(fname,"w");
		fprintf(f,"P5\n%i %i\n255\n",size,size);
		for(int i=0;i<size*size;i++)
			putc((char)(substrate[i].a*255),f);
		fclose(f);

		sprintf(fname,"b%i.ppm",frame);
		f=fopen(fname,"w");
		fprintf(f,"P5\n%i %i\n255\n",size,size);
		for(int i=0;i<size*size;i++)
			putc((char)(substrate[i].b*255),f);
		fclose(f);

	}

}

/*
might wanna optimize diffusion
line by line gaussian?
	this might be good
prolly the best we could do is do multiple iterations in one go

we could also split it among 4 threads?
each covers a quarter
will need 1. threading, 2. job assignments, 3. mutexes
*/

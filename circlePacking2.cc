/*************************************************************************
	> File Name: project3.cc
	> Author: LPQ 
	> Mail: 
	> Created Time: 2016年05月31日 星期二 12时50分48秒
 ************************************************************************/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<float.h>
#include<unistd.h>

#define PRECISION 1e-5
#define LEARNING_RATE 0.1
#define R_RATIO 0.001
#define STEPS 50
#define MAX_FIND 10000
#define MAX_BLOCK 100

typedef struct Circle
{
	float x;
	float y; 
	float radius;
} circle;

float globRadius = 0.24;
float tmpPoint[MAX_FIND], distance[MAX_FIND];

circle blocks[MAX_BLOCK];
int blockNum;

circle getCenterPoint(circle cir[], int n);
void  circleInit(circle **cir, int n);
void getCircleInstance(circle cir[], int circleNum);
void gradientDescent(circle cir[], int n);
void mainSolve(circle userCircle[], int circleNum);
float getCenterDist(circle A, circle B);
float getEnergy(circle cir[], int n);
float sumSquareR(circle cir[], int circleNum);

int main(int argc, char *argv[])
{
    int circleNum;
    clock_t start, end;
    circle *userCircle;
    printf("Please enter the circle numbers\n");

	FILE *fp;
    while (scanf("%d", &circleNum) == 1)
    {
        if (circleNum <= 0)
        {
            printf("the num should be larger than zero!\n");
            continue;
        }

		printf("Please enter the block numbers:\n");
		if (scanf("%d", &blockNum) != 1) {
			printf("Cannot recognize blockNum!\n");
			continue;
		}

		for (int i = 0;i < blockNum;++i) {
			scanf("%f %f", &blocks[i].x, &blocks[i].y);
			blocks[i].radius = 0;
		}

		for (int i = 0;i < blockNum;++i) {
			printf("%f %f\n", blocks[i].x, blocks[i].y);
		}

        start = clock();
        userCircle = (circle*)malloc(circleNum * sizeof(circle));
	    memset(userCircle, 0, circleNum * sizeof(circle));
        //Init the first circle
        //userCircle[0].x = userCircle[0].y = 0.0f;
        //userCircle[0].radius = 1.0f;
        
        mainSolve(userCircle, circleNum);
        end = clock() - start;

		fp = fopen("output.txt", "w");
        for (int i = 0; i < circleNum; i++) {
            printf("circle %d : x: %f y: %f radius: %f energy: %f\n", i, userCircle[i].x, userCircle[i].y, userCircle[i].radius, getEnergy(userCircle, i));
			if (fp)
				fprintf(fp, "%f %f %f\n", userCircle[i].x, userCircle[i].y, userCircle[i].radius);

		}
        printf("sum r^2 : %f time: %f ", sumSquareR(userCircle, circleNum), (float)(end) / CLOCKS_PER_SEC);
        printf("continue or input q to quit\n");
		
		fclose(fp);
    }

	return 0;
}

void circleInit(circle *userCircle, int circleNum)
{
	userCircle = (circle*)malloc(circleNum * sizeof(circle));
	memset(userCircle, 0, circleNum * sizeof(circle));
	
    userCircle[0].x = userCircle[0].y = 0.0f;
    userCircle[0].radius = 1.0f;
 
    return;
}

float sumSquareR(circle userCircle[], int circleNum)
{
    float sum; 
    for(int i = 0; i < circleNum; i++)
        sum += userCircle[i].radius * userCircle[i].radius;
    return sum;
}

float getCenterDist(circle A, circle B)
{
	return sqrt(pow((A.x - B.x), 2) + pow((A.y - B.y), 2));
}

float getEnergy(circle userCircle[], int n)
{
	float energy = 0.0f;
    circle A, B;
    //Energy between each pair of circle
	for (int i = 0; i < n; i++)
	{
		A = userCircle[i];
		for (int j = i + 1; j < n; j++)
		{
            B = userCircle[j];
			if (getCenterDist(A, B) < A.radius + B.radius)
				energy += pow(A.radius + B.radius - getCenterDist(A, B), 2);
		}

		for (int j = 0;j < blockNum;j++) {
			if (getCenterDist(A, blocks[j]) < A.radius)
				energy += pow(A.radius - getCenterDist(A, blocks[j]), 2);
		}
	}
    //Energy to the margins
    for(int i = 0; i < n; i++)
    {
        A = userCircle[i];
		if (A.x + A.radius > 1)
			energy += pow((A.x + A.radius - 1), 2);
		if (A.x - A.radius < -1)
			energy += pow((A.x - A.radius + 1), 2);
		if (A.y + A.radius > 1)
			energy += pow((A.y + A.radius - 1), 2);
		if (A.y - A.radius < -1)
			energy += pow((A.y - A.radius + 1), 2);
    }


    return energy;
}

circle getCenterPoint(circle userCircle[], int circleNum)
{
    int isFind;
    float x , y;
    circle tmpCircle;
    
    isFind = 0;
	if (circleNum == 0) {
        x = 2 * (float)rand() / RAND_MAX - 1;
        y = 2 * (float)rand() / RAND_MAX - 1;
		tmpCircle.x = x;
		tmpCircle.y = y;
		return tmpCircle;
	}

    while(isFind != circleNum)
    {
        isFind = 0;
        x = 2 * (float)rand() / RAND_MAX - 1;
        y = 2 * (float)rand() / RAND_MAX - 1;
    
        for(int i = 0; i < circleNum; i++)
        {
            if(sqrt(pow((x - userCircle[i].x), 2) + pow((y - userCircle[i].y), 2)) > userCircle[i].radius)
            {
                isFind++;
            }
        }
    }

    tmpCircle.x = x;
    tmpCircle.y = y;
    //printf("lpq x %f y %f\n",tmpCircle.x,tmpCircle.y);
    return tmpCircle;

    /*
    tmpCircle.x = tmpCircle.y = 0.0f;
    while(tmpCircle.x == 0.0f || tmpCircle.y == 0.0f)
    {
        memset(tmpPoint, 0.0f, MAX_FIND * sizeof(float));
        memset(distance, 0.0f, MAX_FIND * sizeof(float));
        for(int i = 0; i < MAX_FIND; i++)
        {
            tmpPoint[i] = 2 * (float)rand() / RAND_MAX - 1;
        }
    
        for(int i = 0; i < MAX_FIND - 1; i++)
        {
            x = tmpPoint[i];
            y = tmpPoint[i+1];

            for(int j = 0; j < circleNum; j++)
            distance[i] = sqrt(pow((x - userCircle[j].x), 2) + pow((y - userCircle[j].y), 2)) > userCircle[j].radius ? 1 : 0;
        }

        for(int i = 0; i < MAX_FIND - 1; i++)
        {
            if(distance[i] == 1)
            {
                tmpCircle.x = tmpPoint[i];
                tmpCircle.y = tmpPoint[i+1];
                printf("lpq x %f y %f\n",tmpCircle.x,tmpCircle.y);
                return tmpCircle;
            }
        }

        tmpCircle.x = tmpCircle.y = 0.0f;
    }*/

}

void getCircleInstance(circle userCircle[], int circleNum)
{
    circle MinEnergyCircle, A, B;
    float MinEnergy = FLT_MAX;
    float energy;

    //globRadius = exp(-globRadius);
    //globRadius = globRadius * 0.5f;
    for(int i = 0; i < STEPS; i++)
    {
        energy = 0.0f;
        MinEnergyCircle = getCenterPoint(userCircle, circleNum);
        MinEnergyCircle.radius = globRadius;
        //userCircle[circleNum] = MinEnergyCircle;
        
        A = MinEnergyCircle;
	    for (int j = 0; j < circleNum; j++)
	    {
            B = userCircle[j];
	        if (getCenterDist(A, B) < A.radius + B.radius)
	            energy += pow(A.radius + B.radius - getCenterDist(A, B), 2);
        }
        if (A.x + A.radius > 1)
        	energy += pow((A.x + A.radius - 1), 2);
        if (A.x - A.radius < -1)
        	energy += pow((A.x - A.radius + 1), 2);
        if (A.y + A.radius > 1)
        	energy += pow((A.y + A.radius - 1), 2);
        if (A.y - A.radius < -1)
        	energy += pow((A.y - A.radius + 1), 2);
        if(energy < MinEnergy)
        {
            MinEnergy = energy;
            userCircle[circleNum] = MinEnergyCircle;
        }
    }
    
    return;
}
void gradientDescent(circle userCircle[], int n)
{
    circle A, B;
    A = userCircle[n - 1];
	float x1 = A.x;
	float y1 = A.y;
    float r1 = A.radius;
    printf("gradientDescent\n");
    printf("x1 %f y1 %f radius %f energy: %f\n",x1,y1,r1, getEnergy(userCircle, n));
    float xGrad = 0;
	float yGrad = 0;
    float rGrad = 0;

	for (int i = 0;i < blockNum;++i)
		printf("block%d: %f %f %f\n", i, blocks[i].x, blocks[i].y, blocks[i].radius);

    do {
		xGrad = 0;
		yGrad = 0;
		rGrad = 0;
		A = userCircle[n - 1];
        for (int i = 0; i < n-1; i++)
        {
            B = userCircle[i];
            if (B.radius + r1 - getCenterDist(B, A) > 0)
            {
                xGrad -= 2 * (B.radius + r1 - getCenterDist(B, A)) /
                                getCenterDist(B, A) *
                                (x1 - B.x);
                yGrad -= 2 * (B.radius + r1 - getCenterDist(B, A)) /
                                getCenterDist(B, A) *
                                (y1 - B.y);
                rGrad += 2 * (B.radius + r1 - getCenterDist(B, A));
            }
        }

		for (int j = 0;j < blockNum;j++) {
			B = blocks[j];
			if (r1 > getCenterDist(B, A)) {
                xGrad -= 2 * (r1 - getCenterDist(B, A)) /
                                getCenterDist(B, A) *
                                (x1 - B.x);
                yGrad -= 2 * (r1 - getCenterDist(B, A)) /
                                getCenterDist(B, A) *
                                (y1 - B.y);
                rGrad += 2 * (r1 - getCenterDist(B, A));
			}
		}

        if (r1 + x1 > 1) {
            rGrad += 2 * (x1 + r1 - 1);
            xGrad += 2 * (x1 + r1 - 1);
        }
            
        if (x1 - r1 < -1) {
            rGrad += 2 * (r1 - x1 - 1);
            xGrad -= 2 * (r1 - x1 - 1);
        }
            
        if (y1 + r1 > 1) {
            rGrad += 2 * (y1 + r1 - 1);
            yGrad += 2 * (y1 + r1 - 1);
        }
            
        if (y1 - r1 < -1) {
            rGrad += 2 * (r1 - y1 - 1);
            yGrad -= 2 * (r1 - y1 - 1);
        }

            
		rGrad -= R_RATIO * 2 * r1;

        x1 = x1 - LEARNING_RATE*xGrad;
        y1 = y1 - LEARNING_RATE*yGrad;
        r1 = r1 - LEARNING_RATE*rGrad;
		userCircle[n - 1].x = x1;
		userCircle[n - 1].y = y1;
		userCircle[n - 1].radius = r1;
        //printf("%f %f %f %f\n",x1,y1,r1,getEnergy(userCircle, n));
		//printf("%f %f %f\n", xGrad, yGrad, rGrad);
		//sleep(1000);
    } while (fabs(xGrad) > PRECISION || fabs(yGrad) > PRECISION || fabs(rGrad) > PRECISION);

	userCircle[n - 1].x = x1;
	userCircle[n - 1].y = y1;
	userCircle[n - 1].radius = r1;
    return;
}

void mainSolve(circle userCircle[], int circleNum)
{
    
    for(int i = 0; i < circleNum; i++)
    {
        getCircleInstance(userCircle, i);
        gradientDescent(userCircle, i + 1);
    }
    return;
}


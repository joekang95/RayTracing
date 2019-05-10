# RayTracing

### 0. Running the Program

**For Mac:**

```
> cd hw3
> make
> ./hw3 XXX.scene
```

**For Windows (VS 2017):**

```
> Run the .sln
> Select Project > hw3 Properties > Debugging > Command Arguments
> F5 to Run
```

### 1. Scene Import 

This program supports:

- Single Scene File Input

  ```
  e.g. test2.scene 
  ```

### 2. Input Parameters

	1: Scene File						<XXX.scene>			(required)
	2: Output Image Name				<XXX.jpg>			(optional, but previous parameter need to exist)
	3: Antialiasing 					<Y/y>				(optional, but previous parameter need to exist)
	4: Soft Shadow 						<Y/y>				(optional, but previous parameter need to exist)
	5: Recusrsive Reflection (Sphere) 	<Y/y>				(optional, but previous parameter need to exist)
	6: Animation Output 			  	<Y/y>				(optional, but previous parameter need to exist)
	7: Motion Blur Output 				<Y/y>				(optional, but previous parameter need to exist)


	Example Inputs
	./hw3 test2.scene 			 	  		(No Extra Features)
	./hw3 test2.scene output.jpg 	  		(Output Image)
	./hw3 test2.scene output.jpg y  		(Antialising)
	./hw3 test2.scene output.jpg y y 		(Antialising and Soft Shadow)
	./hw3 test2.scene output.jpg y n y		(Antialising and Recursive Reflection)
	./hw3 test2.scene output.jpg n y n y	(Soft Shadow and Animation)
	./hw3 test2.scene output.jpg y y n n y	(Antialising, Soft Shadow and Motion Blur)

### 3. Keyboard Function

Keyboard ```esc``` -> Exit the Program

### 4. Features

- Ray tracing triangles
- Ray tracing sphere
- Triangle Phong Shading 
- Sphere Phong Shading
- Shadows rays
- Still images 
- Recursive reflection
- Antialiasing
- Soft Shadows
- Animation (Still)
- Motion Blur (Still)

### 5. Still Images

- 001.jpg: test1.scene 
- 002.jpg: test2.scene + antiailasing + soft shadow
- 003.jpg: test2.scene
- 004.jpg: test2.scene + antiailasing + soft shadow + reflection
- 005.jpg: spheres.scene
- 006.jpg: spheres.scene + antiailasing + soft shadow
- 007.jpg: table.scene
- 008.jpg: table.scene + antiailasing + soft shadow
- 009.jpg: SIGGRAPH.scene
- 010.jpg: SIGGRAPH.scene + antiailasing
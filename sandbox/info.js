
rooms.info = function() {

    description = `<b>About the project.</b>
               `;

    code = {

        'Motivation':`

S.html(\`

<b>Motivation</b>
<blockquote>
For my final project I wanted to create a sandbox world 
resembling the beginnings of minecraft for nostalgic value. 
It was one of those games that really affected me because 
it was a place to be creative and learn to work around limitations. 
I enjoyed creative minecraft, so that's the functionality 
I wish to create.
</blockquote>
\`);
`,

        'Application':`

S.html(\`

<b>Application</b>
<blockquote>
To get to the final result, I added a function to texture many different 
types of blocks. I'm procedurally generating the chunks using noise,
and then, to be able to render a seemingly large space while moving
around, I'm only rendering the blocks that can be seen, 
in addition to enabling face culling.
<br>
There are two types of blocks, 
cube blocks and plant blocks. The plant block made up of 2 plane meshes
cross-sectioned.
<br>
<br>
I struggled a lot with the interaction part of the project. The first approach
I took was to render the scene twice, getting the pixel from the screen and reading
values that would describe the block's x, y, z position. However, I could not figure
out how to draw the scene twice and get the pixel only from the first pass.
<br>
Working around these limitations, I decided to pivot to using keyboard keys instead.
Wtih the arrow keys selecting in the x,z plane, and shift/space selecting in the y axis.
This worked, but I definitely felt like it was less intuitive. This is definietely
something I would work in the future.
</blockquote>
\`);
setDescription(description + \`
   <p>
   <table><tr>
   <td><img src=building.png width=350>
   <td width=10>
   <td width=300><big><big><big>
       <b>Building blocks</b>
       <p>
       Demonstration of placed blocks in sandbox.
   </tr></table>

\`);
`,

        'Future':`

S.html(\`

<b>Future</b>
<blockquote>
In the future I would like to add rendering aspects like shadows and water.
I'd also like to add greater variation in the procedural generation, maybe 
ith environment details like trees, and the possibility for cliffs and plains.
<br>
<br>
Again, there were aspects of the interaction that I would expand on and improve
for the future, like adding a visual indication of which block the user has selected,
and make highlighting blocks more intuitive by fixing the x and z axis highlighting
 according to camera pos, rather than world pos. Overall, I had
a lot of fun and learned a lot from this project.

</blockquote>
\`);
`,

    };

}


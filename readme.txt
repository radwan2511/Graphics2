radwan ganem 322509951
moslem asaad 315646802
here we will document our changes/steps we did in the engine:

steps:
same engine from assignment 1 with all the fixes from the privious tirguls
program implemented in main.cpp
implemntaion:
1. get user input (0-5) for desired scene text (list of options is printed in the command line) all scene files are in assignment2 folder.
2. parse the scene text file that the user chose and save its values to the relivant variables from the list below. also print the text file lines.
3. get user input (0-1) option to run with or without anti aliasing (also must have fourth variable in "e" line in the text file as 1.0 to run with anti aliasing).
4. implemented the algorithm like in lecture 3 page 9 and 10 full explanation and documentation will be found in the comments in the code
5. load the image data to the engine and display the result 

Details:

classes implemented:
class ObjectValues {
public:
	int objectIndex;
	float objectColors[4];
	float objectValues[4];
	int objectType;
};

Variables:
# define M_PI           3.14159265 // https://stackoverflow.com/questions/1727881/how-to-use-the-pi-constant-in-c
std::vector<float> camera; // [0] - x , [1] - y , [2] - z , [3] - anti aliasing
std::vector<std::vector<float>> lightDirctions; // [0] - x , [1] - y , [2] - z , [3] - 0.0 for directional light and 1.0 for spotlight.
std::vector<float> ambientLight; // [0] - R , [1] - G , [2] - B , [3] - A
std::vector<std::vector<float>> positions; // [0] - x , [1] - y , [2] - z , [3] - w
std::vector<std::vector<float>> lightIntensities; // [0] - R , [1] - G , [2] - B , [3] - A
//          values, type(object=0,reflictive=1,tranparent=2).     x,y,z center position, r- radius always > 0.   coefficients of the plane equation when 'd' gets is a non - positive value
std::vector<std::pair<std::vector<float>, int>> spheresAndPlanes; // (x,y,z,r) for sphers ,                         (a,b,c,d) for planes
std::vector<std::vector<float>> colors; // [0] - R, [1] - G, [2] - B, [3] - A represents the shininess value
const int DISPLAY_WIDTH = 800;
const int DISPLAY_HEIGHT = 800;


Functions Implemented: 
note: full explation of these function and there implementation are documented in the code.

// spliting each line to a pair of <letter, vactor of the remaining floats>
std::pair<std::string, std::vector<float>> Split(std::string line, char delimiter);

// tirgul 4 pages 4 and 5
/*
Ray-Sphere Intersection I
על מנת לחשב חיתוך של קרן עם כדור, נבצע את החישובים הבאים:
1. נייצג נקודה 𝑃 על הקרן באמצעות משוואת הקרן, ובעזרת 𝑡.
2. נציב במשוואת הכדור את משוואת הנקודה 𝑃.
3. נפתור את המשוואה, וניעזר בנוסחת השורשים.
4. ניקח את הפתרונות עבורם 𝑡>0 (כיוון שעבור ערכים אלו הנקודות נמצאות עם כיוון הקרן).
5. נשים לב כי ייתכנו 3 מקרים:
5.1. לא נקבל אף פיתרון ואז אין חיתוך.
5.2. נקבל פיתרון יחיד ואז יש חיתוך יחיד.
5.3. נקבל 2 פיתרונות ואז יש 2 נקודות חיתוך.
*/
// sphere -> O
// camera -> P_0
// rayCasted -> V
float IntersectSphere(std::vector<float> sphere, std::vector<float> rayCasted, std::vector<float> camera);

// tirgul 4 pages 10 and 11
/*
Ray-Plane Intersection
על מנת לחשב חיתוך של קרן עם מישור, נבצע את החישובים הבאים:
1. נייצג נקודה 𝑃 על הקרן באמצעות משוואת הקרן, ובעזרת 𝑡.
2. נציב במשוואת המישור את משוואת הנקודה 𝑃.
3. נפתור את המשוואה.
4. נשים לב כי ייתכנו 2 מקרים:
4.1. לא נקבל אף פיתרון ואז אין חיתוך.
4.2. נקבל פיתרון יחיד ואז יש חיתוך יחיד.
*/
// Plane -> N
// camera -> (Q_0 - P_0)
// rayCasted -> V
// t = N * (Q_0 - P_0) / N*V
float IntersectPlane(std::vector<float> plane, std::vector<float> rayCasted, std::vector<float> camera);


// see https://computergraphics.stackexchange.com/questions/8479/how-to-calculate-ray
// also tirgul 4 page 21 and page 22
/*
בכדי לחשב את התמונה על המישור שלנו, נבצע את הצעדים הבאים:
1. נגדיר את גודל המישור שלנו, ואת כמות הפיקסלים שלו באורך וברוחב.
2. נחשב האורך והרוחב של פיקסל בודד במישור.
3. נירה קרן אל המרכז של פיקסל ראשוני ממנו אנחנו רוצים להתחיל, ונחשב באיזה צורות הוא פוגע כדי לבחור לו צבע.
4. נזוז למרכזים של כל שאר הפיקסלים במסך (באמצעות קפיצה באורך והרוחב שחישבנו עבור פיקסל בודד), ונירה מהם קרן כדי לבדוק באילו צורות פגענו ואיזה צבע לבחור להם. 	
*/
// P = P_c + ( x - | R_x/2 | ) R V_right - ( y - | R_y/2 | ) R V_up
// note from assignment pdf: 
// The screen is located on z=0 plane. The right up corner of the screen is located at (1,1,0) and the
// left bottom corner of the screen is located at(-1, -1, 0).All in the scene coordinates.	
std::vector<float> ConstructRayThroughPixel(int y, int x);

// page 19 in lecture 3 algorithm
std::pair<float[4], ObjectValues> FindIntersection(std::vector<float> ray, int objectIndexToPass, std::vector<float> camera);

// Get Color of Object / Reflective Object / Transparent Object
std::vector<float> GetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int times, std::vector<float> camera);

// Calculates Diffuse Reflection, Specular Reflection And Shadow Term and update the color
std::vector<float> ObjectGetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, std::vector<float> camera);

// if it is sphere simply return the color
// if it is plane calculate color same as tirgul 5 page 4
std::vector<float> GetObjectColor(ObjectValues obj, float* hit);

// tirgul 4 page 24
// calculate diffuse reflection
// cosine law 
// N * L = |N||L|cos a
// N^ * L^ = cos a
// I_D = K_D ( N^ * L^ ) I_L
/*
על מנת לחשב את ההחזרה הדיפיוזית, נבצע את החישובים הבאים:
1. נחשב את המכפלה הסקאלרית של וקטור הנורמל למישור עם הוקטור מהמישור למקור האור.
2. נכפיל את התוצאה במקדם הדיפיוזי 𝐾_𝐷 ובעוצמת מקור האור 𝐼_𝐿.
3. את התוצאה נסמן כ- 𝐼_𝐷.
*/
std::vector<float> DiffuseReflection(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int i, std::vector<float> camera);

// calculate Specular Reflection
// tirgul 4 pages 26 + 27
// phong model
// cos(a)^n
// V^ * R^ = max( 0 , V^ * R^ )
// I_S = K_S ( V^ * R^ )^n I_L
/*
על מנת לחשב את החזר הראי, נבצע את החישובים הבאים:
1. נחשבת את המכפלה הסקאלרית בין הוקטור מהמישור למקור האור עם הוקטור מהמישור אל הצופה (או במקרה שלנו אל המצלמה), וניקח את המקסימום בין ערך זה ל- 0.
2. את התוצאה נכפיל בחזקת מקדם הברק 𝑛.
3. נכפיל את התוצאה במקדם הראי 𝐾_𝑆 ובעוצמת מקור האור 𝐼_𝐿.
4. את התוצאה נסמן כ- 𝐼_𝑆.
*/
std::vector<float> SpecularReflection(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int i, std::vector<float> camera);

// calculate Shadow term
// tirgul 4 pages 31 and 32
/*
מקדם הצל אומר לנו איזה מקורות אור חסומים עבור הצורה שלנו (כלומר קרן האור לא מצליחה לפגוע בהם).
לכן בכדי לחשב את הצל על הצורה שלנו, כל סכום של ההחזרה הדיפיוזית והחזר הראי, נכפיל במקדם הצל 𝑆_𝑖.
מקדם הצל יהיה שווה ל- 0 אם מקור האור חסום עבור הצורה שלנו.
אחרת, מקדם הצל יהיה שווה ל- 1.
*/
std::vector<float> Shadows(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int i, std::vector<float> camera);

/*
Trace secondary ray in direction of mirror reflection
Evaluate radiance along secondary ray and include it into illumination model
במידה ויש לנו חומר מראתי, נחשב את הצבע שלו על ידי כך שנוסיף לסכום את החישוב הרקורסיבי של פגיעת קרן האור בגופים נוספים לאחר שהיא נשברה והמשיכה מחוץ לצורה שלנו.
*/
std::vector<float> ReflectiveGetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int times, std::vector<float> camera);

/*
במידה ויש לנו חומר שקוף, נחשב את הצבע שלו על ידי כך שנוסיף לסכום את החישוב הרקורסיבי של פגיעת קרן האור בגופים נוספים לאחר שהיא נשברה ועברה דרך הצורה שלנו.
*/
// https://personal.math.ubc.ca/~cass/courses/m309-01a/chu/Fundamentals/snell.htm
// https://www.physicsclassroom.com/class/refrn/Lesson-2/Snell-s-Law
// from assignment pdf: Use Snell's law with refractive index of 1.5 to determine the direction of the light rays intersect
// with the sphere(air have refractive index of 1)
// snell's law
// n_r * sin a_r = n_i * sin a_i
// T = ( ( n_i / n_r ) * cos a_i - cos a_r ) * N - ( n_i / n_r ) * L )
std::vector<float> TranparentGetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int times, std::vector<float> camera);


// see https://computergraphics.stackexchange.com/questions/8479/how-to-calculate-ray
// also see explanation in the ConstructRayThroughPixel explanation
// shoting 4 rays: top left, top right, buttom left, buttom right	
std::vector<std::vector<float>> ConstructRayThroughPixelBonus(int y, int x);


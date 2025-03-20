from PIL import Image, ImageDraw
def draw(n, m):
    path = "/Users/daniil/Documents/lab/ray_tracing/res.txt"
    img = Image.new(mode = "RGB", size = (n, m))
    draw = ImageDraw.Draw(img, mode = "RGB")
    f1 = open(path, 'r')
    b = (f1.read().split())
    f1.close()
    x_r = 10**(-20)
    for i in range(n):
        for j in range(m):
            cur = float(b[3 * (i * m + j)])
            cur = max(cur, float(b[3 * (i * m + j) + 1]))
            cur = max(cur, float(b[3 * (i * m + j) + 2]))
            if (x_r < cur):
                x_r = cur
    print (x_r)
    for i in range(n):
        for j in range(m):
            rq = 1
            a_r = int(float(b[3 * (i * m + j)]) * rq * 255.0)
            a_g = int(float(b[3 * (i * m + j) + 1]) * rq * 255.0)
            a_b = int(float(b[3 * (i * m + j) + 2]) * rq * 255.0)
            draw.point((i, m - 1 - j), fill = (min(a_r, 2550), min(a_g, 2550), min(a_b, 2550)))
    img.show()
path = "/Users/daniil/Documents/lab/ray_tracing/conf.txt"
f0 = open(path, 'r')
b = (f0.read().split())
f0.close()
n = int(b[0])
m = int(b[1])
print(n, m)
draw(n, m)

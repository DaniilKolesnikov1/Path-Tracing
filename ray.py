from PIL import Image, ImageDraw
def draw(n, m):
    path = "/Users/daniil/Documents/lab/ray_tracing/res.txt"
    img = Image.new(mode = "RGB", size = (n, m))
    draw = ImageDraw.Draw(img, mode = "RGB")
    f1 = open(path, 'r')
    b = (f1.read().split())
    f1.close()
    x_r = 10**(-20)
    x_q = 1
    for i in range(n):
        for j in range(m):
            cur = float(b[3 * (i * m + j)])
            cur = max(cur, float(b[3 * (i * m + j) + 1]))
            cur = max(cur, float(b[3 * (i * m + j) + 2]))
            cur = cur
            x_q = min(x_q, float(b[3 * (i * m + j)]))
            if (x_r < cur):
                x_r = cur
    print (x_r, x_q)
    for i in range(n):
        for j in range(m):
            rq = 1
            a_r = int(float(b[3 * (i * m + j)]) * rq * 255.0 / x_r)
            a_g = int(float(b[3 * (i * m + j) + 1]) * rq * 255.0 / x_r)
            a_b = int(float(b[3 * (i * m + j) + 2]) * rq * 255.0 / x_r)
            draw.point((i, m - 1 - j), fill = (min(a_r, 255), min(a_g, 255), min(a_b, 255)))
    img.show()
path = "/Users/daniil/Documents/lab/ray_tracing/conf.txt"
f0 = open(path, 'r')
b = (f0.read().split())
f0.close()
n = int(b[0])
m = int(b[1])
draw(n, m)

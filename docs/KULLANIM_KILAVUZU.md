# CoreForge — Kullanım Kılavuzu

Modüler reaktör (SMR/MMR) tasarımını destekleyen **2D/3D çok-gruplu
nötronik analiz ve simülasyon kod sistemi**. Bu kılavuz kurulumdan tüm
sekmelerin kullanımına kadar adım adım rehberdir (Teknofest NDK
dokümantasyon şartı kapsamında).

## 1 · Kurulum

Gereksinimler: Linux/WSL2, bir Fortran derleyicisi (Intel `ifx` veya
`gfortran`), Python ≥ 3.9.

```bash
# 1) Fortran motorunu derle
ifx -O3 -qopenmp solver/coreforge.f90 -o solver/coreforge
#    (alternatif: gfortran -O3 -fopenmp ...)
#    Bu adım atlanabilir: motor yoksa Python sürücüsü ilk kullanımda
#    ifx/gfortran ile KENDİSİ derler (günlük: solver/build.log).

# 2) Python bağımlılıkları
pip install -r requirements.txt

# 3) DOĞRULAMA — kuruluma her zaman bununla başlayın
python3 verify.py            # 24 kontrol, ~10 s, hepsi PASS olmalı
#    Fortran derleyicisi olmayan ortamda: python3 verify.py --no-engine

# 4) Arayüzü başlat
streamlit run app.py         # tarayıcı: http://localhost:8501
```

> Canlı demo / bulut dağıtımı için `docs/YAYINLAMA.md`, bağımsız
> doğrulama turu için `docs/BAGIMSIZ_DOGRULAMA.md`.

`verify.py` çıktısındaki her satır, DOGRULAMA_RAPORU.md'de açıklanan bir
benchmark/analitik kontroldür. Tek bir FAIL bile kurulumun sağlıksız
olduğunu gösterir.

## 2 · Genel iş akışı

```
Benchmark yükle VEYA sıfırdan kur → (yakıt tasarla) → kor deseni →
2D/3D seç → ÇÖZ → sonuçları incele → (çevrim/transient/araçlar) → rapor
```

Sol kenar çubuğu her zaman görünürdür: çözücü modeli (2-D/3-D), mesh
inceltme, sınır koşulları, gelişmiş ayarlar ve **proje kaydet/yükle**.

## 3 · Sekmeler

### 📚 Benchmarks
Tek tıkla yüklenen 7 hazır kor (SMR, IAEA-3D, IAEA-2D, C5G7, Designer
PWR, mini-kor, k∞). Kartlardaki 🧊 3-D / ▦ 2-D etiketi çözücü modelini
gösterir. Yüklenen her şey düzenlenebilir bir başlangıç noktasıdır.
Alttaki doğrulama tablosu bu motorun ölçülmüş sonuçlarıdır.

### 🧬 Fuel designer
Fiziksel girdiden (U-235 zenginliği, çözünür boron ppm, moderatör
yoğunluğu, pin yarıçapı) 2-gruplu makroskopik kesit üretir; k∞ canlı
gösterilir. **Add to materials** ile kor malzemesi olur.
Alt bölümdeki **🔁 Equivalent fuel**, yalnızca kesitle tanımlı bir
yakıta (ör. benchmark yakıtı) eşdeğer kompozisyon oturtur — böylece o
malzeme yakılabilir (burnup) hale gelir; kalan farklar tablo halinde
dürüstçe raporlanır.

### 🧪 Materials
Grup sayısı (1–8), malzeme sayısı (≤8), D/Σa/νΣf/χ tablosu ve her
malzeme için tam saçılma matrisi (up-scatter serbest). Benchmark
sabitleri ya da lattice-code çıktısı doğrudan buraya yazılabilir.

### 🧱 Core builder
- **2-D**: pitch, satır/sütun, ızgaraya malzeme numarası; canlı renkli
  önizleme. Kenar çubuğundaki anahtarla kor 3-D'ye "extrude" edilir.
- **3-D**: eksenel zon yöneticisi — zon ekle/sil/taşı/kopyala, katman
  sayısı, dz, **▩ fill** ve **↷ copy** araçları; zon seçicide fiziksel
  z-aralıkları görünür. Her zonun kendi radyal haritası vardır
  (kısmi batmış çubuklar, eksenel reflektörler böyle kurulur).

### ⚡ Solve & results
**Solve k_eff** → metrikler (k_eff, ρ, referans farkı, F_q, 3-D'de F_z),
radyal harita ve güç/akı haritaları, 3-D'de **eksenel güç profili** ve
katman gezgini, assembly güç haritası, akı traversi, nötron dengesi.
**🔌 Operating point** paneli, verilen kor termal gücünden (MW) mutlak
birimleri üretir: W/cm³ güç yoğunluğu, n/cm²·s akı, assembly başına MW.
İndirmeler: input.txt, flux.csv, power.csv ve **📄 tek dosyalık HTML
rapor** (QA/izlenebilirlik bölümü dâhil). Oturum içi koşu geçmişi en
alttadır.

**🔒 QA bütünlük kilidi:** Bir benchmark yüklendikten sonra geometri,
malzeme veya sınır koşulu değiştirilirse "Δ vs reference" otomatik
olarak *invalidated* olur ve uyarı gösterilir — yayınlanmış referans
artık o problemi tanımlamaz çünkü. Mesh inceltme (div/divz) referansı
BOZMAZ; yakınsama çalışmaları serbesttir. Preset'i yeniden yüklemek
kilidi sıfırlar.

### 🌡 Thermal-hydraulics
Nötronik çözümden beslenen kapalı-kanal analiz: ortalama ve **sıcak
kanal** için soğutucu ısınması, kılıf ve **yakıt merkez sıcaklığı**
(Dittus-Boelter + iletim dirençleri), sürtünme **basınç düşümü**, doyma
marjı. Radyal tepe (F_ΔH) güç haritasından, eksenel şekil 3-D profilden
gelir. Tek-fazlı model; kaynama başlangıcı marj ile işaretlenir.

### 🔥 Burnup
Blok-blok tükenim (U→Pu zinciri, denge Xe/Sm, toplu FP), opsiyonel
**kritik boron letdown** ve **otomatik çevrim sonu** (boron 0 ppm'e
düşüp kritiklik korunamayınca çevrim biter → ulaşılabilir yanma +
EFPD). **Çok-çevrim:** çevrim sonunda 🔁 en yanmış blokları taze yakıtla değiştirip (batch fraction) sonraki çevrimi koşabilirsin — denge çevrimine yaklaşımı izle. Çıktılar: k(BU), letdown eğrisi, EOC yanma haritası, güç kayması
haritası, CSV. (Bu sürümde 2-D korlar için; designer/eşdeğer yakıt
gerektirir.)

### ⏱ Transient
Nokta kinetiği (6 gecikmiş grup) ile zamana bağlı güç: **step / rampa /
scram** senaryoları, opsiyonel **Doppler geri beslemesi** (α_D, yakıt
ısı kapasitesi, ısı çekme zaman sabiti). Donmuş-ρ matris-üstel adımıyla
kararlılık sorunu yoktur; asimptotik periyot analitik **inhour**
denklemiyle karşılaştırılır. Not: bozunum ısısı modellenmez.
Alt bölümde **☁️ Xenon transient**: güç değişimi sonrası Xe-135 reaktivitesi — kapanma sonrası **iyot çukuru** (derinlik/zaman) ve toparlanma süresi.

### 🔧 Physics tools
1. **Malzeme takası değeri** — ör. kontrol çubuğu grubunun pcm değeri
2. **Çubuk batırma taraması (3-D)** — IAEA-3D'de çubuk-5 derinliği ↦
   k eğrisi: klasik S-eğrisi
3. **Kritik çubuk pozisyonu (3-D)** — hedef k için batma derinliği (bisection)
4. **Kritik boron araması [ppm]** — designer korlarda k=1 için boron
5. **Genel ΔΣa kritiklik araması** — herhangi bir malzeme/grupta

## 4 · Proje kaydet/yükle

Kenar çubuğu → **💾 Save project**: tüm konfigürasyon (malzemeler,
nüklit vektörleri, zonlar, BC, ayarlar) tek JSON dosyası.
**📂 Load project** ile geri yüklenir. **🆕 New blank core** temiz
başlangıç verir.

## 5 · Motoru arayüzsüz kullanma

```bash
./solver/coreforge input.txt
```
Girdi biçimi README'de belgelenmiştir (anahtar kelimeli, `#` yorum).
Çıktı: stdout'ta KEFF/FXY/denge; `flux.csv`, `power.csv`.

## 6 · Sık karşılaşılan durumlar

| Belirti | Neden / çözüm |
|---|---|
| "Engine not built" | Otomatik derleme de başarısız — `solver/build.log`'a bakın; ifx/gfortran kurulu mu? (Kurulum adım 1) |
| Solve çok yavaş | mesh kaydırıcılarını düşürün; küçük problemde "auto threads" açık kalsın |
| Burnup "designer fuel yok" | 🧬'den yakıt ekleyin ya da eşdeğer yakıt iliştirin |
| 3-D'de burnup kapalı | Bu sürümde tükenim radyal (2-D) modeldedir |
| k_eff NaN / diverged | Gelişmiş ayarlarda ω'yı düşürün (ör. 1.4) |

## 7 · Kapsam ve sınırlar (dürüst beyan)

Difüzyon teorisi (pin-ölçeği transport değildir — bkz. C5G7 demosu),
donmuş 2-grup spektrum, T-H geri beslemesi yok (transient'te toplu yakıt
sıcaklığı modeli hariç), bozunum ısısı yok, tükenim 2-D. Ayrıntı:
TEORI_VE_YONTEM.md §Sınırlar.

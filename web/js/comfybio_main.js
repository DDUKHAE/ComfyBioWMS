import { app } from "../../scripts/app.js";

// Domain color schemes for ComfyBIO nodes
const DOMAIN_COLORS = {
    "ComfyBIO/Core": { color: "#1b4332", bgcolor: "#2d6a4f", groupcolor: "#081c15" },
    "ComfyBIO/Variant": { color: "#3c096c", bgcolor: "#5a189a", groupcolor: "#240046" },
    "ComfyBIO/ATAC": { color: "#7b2cbf", bgcolor: "#9d4edd", groupcolor: "#3c096c" },
    "ComfyBIO/Metagenome": { color: "#004b23", bgcolor: "#007200", groupcolor: "#002800" },
    "ComfyBIO/Assembly": { color: "#006466", bgcolor: "#065a60", groupcolor: "#0b525b" },
    "ComfyBIO/Visualizer": { color: "#0a9396", bgcolor: "#005f73", groupcolor: "#001219" },
    "ComfyBIO/Biopython": { color: "#ca6702", bgcolor: "#bb3e03", groupcolor: "#9b2226" },
    "ComfyBIO/Genomics": { color: "#1d3557", bgcolor: "#457b9d", groupcolor: "#14213d" },
    "ComfyBIO/SingleCell": { color: "#6a040f", bgcolor: "#9d0208", groupcolor: "#370617" },
    "ComfyBIO/Epigenomics": { color: "#3f37c9", bgcolor: "#4361ee", groupcolor: "#3a0ca3" },
    "ComfyBIO/Proteomics": { color: "#7209b7", bgcolor: "#b5179e", groupcolor: "#480ca8" },
    "ComfyBIO/CADD": { color: "#d00000", bgcolor: "#dc2f02", groupcolor: "#6a040f" },
    "ComfyBIO/Microbiome": { color: "#2b9348", bgcolor: "#55a630", groupcolor: "#1b4332" },
};

app.registerExtension({
    name: "ComfyBIO.WMS",
    async beforeRegisterNodeDef(nodeType, nodeData, app) {
        const category = nodeData.category || "";
        
        // Check if this is a ComfyBIO node
        if (category.startsWith("ComfyBIO") || category.includes("bio") || category.includes("Bio")) {
            // Apply category-specific color if matched
            for (const [prefix, theme] of Object.entries(DOMAIN_COLORS)) {
                if (category.startsWith(prefix)) {
                    nodeType.prototype.color = theme.color;
                    nodeType.prototype.bgcolor = theme.bgcolor;
                    nodeType.prototype.groupcolor = theme.groupcolor;
                    break;
                }
            }

            // Customize node UI appearance on creation
            const onNodeCreated = nodeType.prototype.onNodeCreated;
            nodeType.prototype.onNodeCreated = function () {
                const r = onNodeCreated ? onNodeCreated.apply(this, arguments) : undefined;
                this.badge = "🧬 BIO";
                return r;
            };
        }
    },
    async setup() {
        console.log("%c[ComfyBIOWMS]%c High-throughput Bioinformatics WMS Extension Loaded", "color: #2b9348; font-weight: bold;", "color: inherit;");
    }
});
